PRIVATE unsigned char
hc_int_cb_def(int   i,
              int   j,
              int   k,
              int   l,
              vrna_hc_t *hc);


PRIVATE unsigned char
hc_int_cb_def_user(int  i,
                   int  j,
                   int  k,
                   int  l,
                   vrna_hc_t *hc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE unsigned char
wrap_hc_int_cb(int   i,
              int   j,
              int   k,
              int   l,
              vrna_hc_t *hc)
{
  unsigned char pij, pkl;
  unsigned int  u1, u2, u3;

  if (j < k) {
    /* circular RNA */
    u1 = i - 1;
    u2 = k - j;
    u3 = hc->n - l;

    if ((u1) &&
        (hc->up_int[1] < u1))
      return (unsigned char)0;

    if ((u2) &&
        (hc->up_int[j + 1] < u2))
      return (unsigned char)0;

    if ((u3) &&
       (hc->up_int[l + 1] < u3))
      return (unsigned char)0;

  } else {
    /* linear RNA */
    u1 = k - i - 1;
    u2 = j - l - 1;

    if ((u1) &&
        (hc->up_int[i + 1] < u1))
      return (unsigned char)0;

    if ((u2) &&
        (hc->up_int[l + 1] < u2))
      return (unsigned char)0;
  }

  pij = hc->mx[hc->n * i + j];
  pkl = hc->mx[hc->n * k + l];

  if ((pij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) &&
      (pkl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC))
    return (unsigned char)1;

  return (unsigned char)0;
}


PRIVATE unsigned char
hc_int_cb_sn(int   i,
             int   j,
             int   k,
             int   l,
             vrna_hc_t *hc)
{
  unsigned char pij, pkl;
  unsigned int  u1, u2, u3;

  if ((hc->sn[i] != hc->sn[k]) ||
      (hc->sn[l] != hc->sn[j]))
    return (unsigned char)0;

  return (unsigned char)1;
}


PRIVATE unsigned char
wrap_hc_int_cb_window(int   i,
              int   j,
              int   k,
              int   l,
              vrna_hc_t *hc)
{
  unsigned char pij, pkl;
  unsigned int  u1, u2, u3;

  u1 = k - i - 1;
  u2 = j - l - 1;

  if ((u1) &&
      (hc->up_int[i + 1] < u1))
    return (unsigned char)0;

  if ((u2) &&
      (hc->up_int[l + 1] < u2))
    return (unsigned char)0;

  pij = hc->matrix_local[i][j - i];
  pkl = hc->matrix_local[k][l - k];

  if ((pij & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) &&
      (pkl & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC))
    return (unsigned char)1;

  return (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_int_cb_sn(int            i,
                 int            j,
              int   k,
              int   l,
                 vrna_hc_t      *hc)
{
  return (hc_int_cb_sn(i, j, k, l, hc)) ? wrap_hc_int_cb(i, j, k, l, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_int_cb_user(int            i,
                   int            j,
              int   k,
              int   l,
                   vrna_hc_t      *hc)
{
  return (hc->f(i, j, k, l, VRNA_DECOMP_PAIR_IL, hc->data)) ? wrap_hc_int_cb(i, j, k, l, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_int_cb_sn_user(int           i,
                      int           j,
              int   k,
              int   l,
                      vrna_hc_t     *hc)
{
  return (hc->f(i, j, k, l, VRNA_DECOMP_PAIR_IL, hc->data)) ? wrap_hc_int_cb_sn(i, j, k, l, hc) : (unsigned char)0;
}


PRIVATE unsigned char
wrap_hc_int_cb_user_window(int           i,
                          int           j,
              int   k,
              int   l,
                          vrna_hc_t     *hc)
{
  return (hc->f(i, j, k, l, VRNA_DECOMP_PAIR_IL, hc->data)) ? wrap_hc_int_cb_window(i, j, k, l, hc) : (unsigned char)0;
}
