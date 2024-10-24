/*
 *  This file contains hard constraint wrappers required by both, MFE and
 *  partition function version of exterior loop evaluation
 */

PRIVATE unsigned char
wrap_hc_ext_cb(int           i,
              int           j,
              int           k,
              int           l,
              unsigned char d,
              vrna_hc_t     *hc);


PRIVATE unsigned char
wrap_hc_ext_cb_sn(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             vrna_hc_t      *hc);


PRIVATE unsigned char
wrap_hc_ext_cb_window(int            i,
                     int            j,
                     int            k,
                     int            l,
                     unsigned char  d,
                     vrna_hc_t      *hc);


PRIVATE unsigned char
wrap_hc_ext_cb_user(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   vrna_hc_t      *hc);


PRIVATE unsigned char
wrap_hc_ext_cb_sn_user(int           i,
                      int           j,
                      int           k,
                      int           l,
                      unsigned char d,
                      vrna_hc_t     *hc);


PRIVATE unsigned char
wrap_hc_ext_cb_user_window(int           i,
                          int           j,
                          int           k,
                          int           l,
                          unsigned char d,
                          vrna_hc_t     *hc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE unsigned char
wrap_hc_ext_cb(int           i,
              int           j,
              int           k,
              int           l,
              unsigned char d,
              vrna_hc_t     *hc)
{
  unsigned int          n, di, dj;
  unsigned char         eval;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  n     = hc->n;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      if (hc->mx[n * j + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (hc->up_ext[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT:
      if (hc->mx[n * k + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (hc->up_ext[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (hc->mx[n * (j - 1) + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (hc->up_ext[j] == 0)
          eval = (unsigned char)0;

        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;

          if ((di != 0) && (hc->up_ext[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (hc->mx[n * k + i + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (hc->up_ext[i] == 0)
          eval = (unsigned char)0;

        if (j != k) {
          /* otherwise, stem spans from i + 1 to j */
          dj = l - k - 1;

          if ((dj != 0) && (hc->up_ext[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      eval  = (unsigned char)1;
      di    = l - k - 1;
      if ((di != 0) && (hc->up_ext[k + 1] < di))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_STEM:
      if (hc->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (hc->up_ext[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (hc->up_ext[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (hc->up_ext[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (hc->up_ext[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (hc->up_ext[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    case VRNA_DECOMP_EXT_STEM_OUTSIDE:
      if (hc->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)
        eval = (unsigned char)1;

      break;

    default:
      vrna_log_warning("Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
hc_ext_cb_sn(int            i,
             int            j,
             int            k,
             int            l,
             unsigned char  d,
             vrna_hc_t      *hc)
{
  unsigned int          *sn;
  unsigned char         eval;

  sn    = hc->sn;
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
      vrna_log_warning("Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
wrap_hc_ext_cb_window(int            i,
                     int            j,
                     int            k,
                     int            l,
                     unsigned char  d,
                     vrna_hc_t      *hc)
{
  unsigned int          di, dj;
  unsigned char         eval;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      if (hc->matrix_local[l][j - l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (hc->up_ext[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT:
      if (hc->matrix_local[i][k - i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (j != k) {
          /* otherwise, stem spans from i to j */
          dj = l - k - 1;
          if ((dj != 0) && (hc->up_ext[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      if (hc->matrix_local[l][j - 1 - l] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (hc->up_ext[j] == 0)
          eval = (unsigned char)0;

        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;

          if ((di != 0) && (hc->up_ext[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      if (hc->matrix_local[i + 1][k - (i + 1)] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (hc->up_ext[i] == 0)
          eval = (unsigned char)0;

        if (j != k) {
          /* otherwise, stem spans from i + 1 to j */
          dj = l - k - 1;

          if ((dj != 0) && (hc->up_ext[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM:
      if (hc->matrix_local[k][l - k] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (hc->up_ext[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (hc->up_ext[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT_EXT:
      eval  = (unsigned char)1;
      di    = l - k - 1;
      if ((di != 0) && (hc->up_ext[k + 1] < di))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (hc->up_ext[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (hc->up_ext[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (hc->up_ext[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    default:
      vrna_log_warning("Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
wrap_hc_ext_cb_sn(int            i,
                 int            j,
                 int            k,
                 int            l,
                 unsigned char  d,
                 vrna_hc_t      *hc)
{
  return (hc_ext_cb_sn(i, j, k, l, d, hc)) ? wrap_hc_ext_cb(i, j, k, l, d, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_ext_cb_user(int            i,
                   int            j,
                   int            k,
                   int            l,
                   unsigned char  d,
                   vrna_hc_t      *hc)
{
  return (hc->f(i, j, k, l, d, hc->data)) ? wrap_hc_ext_cb(i, j, k, l, d, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_ext_cb_sn_user(int           i,
                      int           j,
                      int           k,
                      int           l,
                      unsigned char d,
                      vrna_hc_t     *hc)
{
  return (hc->f(i, j, k, l, d, hc->data)) ? wrap_hc_ext_cb_sn(i, j, k, l, d, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_ext_cb_user_window(int           i,
                          int           j,
                          int           k,
                          int           l,
                          unsigned char d,
                          vrna_hc_t     *hc)
{
  return (hc->f(i, j, k, l, d, hc->data)) ? wrap_hc_ext_cb_window(i, j, k, l, d, hc) : (unsigned char)0;
}


PRIVATE unsigned char
dispatch_wrap_ext_f(int            i,
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
      hc->eval_ext   = wrap_hc_ext_cb_user_window;
    } else {
      hc->eval_ext   = wrap_hc_ext_cb_window;
    }
  } else {
    if (hc->strands == 1) {
      if (hc->f) {
        hc->eval_ext   = wrap_hc_ext_cb_user;
      } else {
        hc->eval_ext   = wrap_hc_ext_cb;
      }
    } else {
      if (hc->f) {
        hc->eval_ext   = wrap_hc_ext_cb_sn_user;
      } else {
        hc->eval_ext   = wrap_hc_ext_cb_sn;
      }
    }
  }

  /* call eval wrapper function */
  return hc->eval_ext(i, j, k, l, d, hc);
}
