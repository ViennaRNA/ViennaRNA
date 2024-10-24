

PRIVATE unsigned char
wrap_hc_hp_cb_def(int            i,
             int            j,
             vrna_hc_t *hc);


PRIVATE unsigned char
wrap_hc_hp_cb_def_window(int           i,
                    int           j,
                    vrna_hc_t *hc);


PRIVATE unsigned char
wrap_hc_hp_cb_def_user(int           i,
                  int           j,
                  vrna_hc_t *hc);


PRIVATE unsigned char
wrap_hc_hp_cb_def_user_window(int            i,
                         int            j,
                         vrna_hc_t *hc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE unsigned char
wrap_hc_hp_cb(int       i,
              int       j,
              vrna_hc_t *hc)
{
  unsigned int          u;
  unsigned char         eval;

  eval = (char)0;

  if (hc->mx[hc->n * i + j] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
    if (j > i) {
      /* linear case */
      u = j - i - 1;
    } else {
      /* circular case */
      u = hc->n - i + j - 1;
    }

    if (hc->up_hp[i + 1] >= u)
      eval = (unsigned char)1;
  }

  return eval;
}


PRIVATE unsigned char
hc_hp_cb_sn(int       i,
            int       j,
            vrna_hc_t *hc)
{
  /* no strand nicks are allowed in hairpin loops */
  if (hc->sn[i] != hc->sn[j])
    return (unsigned char)0;

  return (unsigned char)1;
}


PRIVATE unsigned char
wrap_hc_hp_cb_window(int        i,
                     int        j,
                     vrna_hc_t  *hc)
{
  unsigned int          u;
  unsigned char         eval;

  eval = (unsigned char)0;

  u = j - i - 1;

  if (hc->matrix_local[i][j - i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
    eval = (unsigned char)1;
    if (hc->up_hp[i + 1] < u)
      eval = (unsigned char)0;
  }

  return eval;
}


PRIVATE unsigned char
wrap_hc_hp_cb_sn(int            i,
                 int            j,
                 vrna_hc_t      *hc)
{
  return (hc_hp_cb_sn(i, j, hc)) ? wrap_hc_hp_cb(i, j, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_hp_cb_user(int            i,
                   int            j,
                   vrna_hc_t      *hc)
{
  return (hc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, hc->data)) ? wrap_hc_hp_cb(i, j, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_hp_cb_sn_user(int           i,
                      int           j,
                      vrna_hc_t     *hc)
{
  return (hc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, hc->data)) ? wrap_hc_hp_cb_sn(i, j, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_hp_cb_user_window(int           i,
                          int           j,
                          vrna_hc_t     *hc)
{
  return (hc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, hc->data)) ? wrap_hc_hp_cb_window(i, j, hc) : (unsigned char)0;
}


PRIVATE unsigned char
dispatch_wrap_hp_f(int        i,
                   int        j,
                   vrna_hc_t  *hc)
{
  /* check whether additional preparations are necessary, e.g. due to constrahps added in the meantime */
  if (hc->state != STATE_CLEAN) {
    vrna_log_error("uninitialized hard constraints! Call vrna_hc_prepare() first!");
    return (unsigned char)0;
  }

  /* assign eval wrapper functions */
  if (hc->type == VRNA_HC_WINDOW) {
    if (hc->f) {
      hc->eval_hp   = wrap_hc_hp_cb_user_window;
    } else {
      hc->eval_hp   = wrap_hc_hp_cb_window;
    }
  } else {
    if (hc->strands == 1) {
      if (hc->f) {
        hc->eval_hp   = wrap_hc_hp_cb_user;
      } else {
        hc->eval_hp   = wrap_hc_hp_cb;
      }
    } else {
      if (hc->f) {
        hc->eval_hp   = wrap_hc_hp_cb_sn_user;
      } else {
        hc->eval_hp   = wrap_hc_hp_cb_sn;
      }
    }
  }

  /* call eval wrapper function */
  return hc->eval_hp(i, j, hc);
}
