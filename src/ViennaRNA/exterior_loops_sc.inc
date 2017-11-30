/*
 *  This file contains soft constraints wrappers for exterior loop evaluation
 */

struct sc_wrapper_f5;

typedef int (sc_f5_reduce_to_ext)(int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_wrapper_f5  *data);


typedef int (sc_f5_reduce_to_stem)(int                  j,
                                   int                  k,
                                   int                  l,
                                   struct sc_wrapper_f5 *data);


typedef int (sc_f5_split_in_ext_stem)(int                   j,
                                      int                   k,
                                      int                   l,
                                      struct sc_wrapper_f5  *data);


typedef int (sc_f5_split_in_ext_stem1)(int                  j,
                                       int                  k,
                                       int                  l,
                                       struct sc_wrapper_f5 *data);


struct sc_wrapper_f5 {
  int                       **up;

  sc_f5_reduce_to_ext       *red_ext;
  sc_f5_reduce_to_stem      *red_stem;
  sc_f5_split_in_ext_stem   *decomp_stem;
  sc_f5_split_in_ext_stem1  *decomp_stem1;

  vrna_callback_sc_energy   *user_cb;
  void                      *user_data;

  /* below attributes are for comparative structure prediction */
  int                       n_seq;
  unsigned int              **a2s;
  int                       ***up_comparative;

  vrna_callback_sc_energy   **user_cb_comparative;
  void                      **user_data_comparative;
};


PRIVATE INLINE void
init_sc_wrapper(vrna_fold_compound_t  *fc,
                struct sc_wrapper_f5  *sc_wrapper);


PRIVATE INLINE void
free_sc_wrapper(struct sc_wrapper_f5 *sc_wrapper);


PRIVATE INLINE int
sc_default_f5_reduce(int                  j,
                     int                  k,
                     int                  l,
                     struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_def_f5_reduce_to_ext(int                  j,
                             int                  k,
                             int                  l,
                             struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_f5_reduce_to_ext(int                  j,
                         int                  k,
                         int                  l,
                         struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_default_f5_reduce_comparative(int                  j,
                                 int                  k,
                                 int                  l,
                                 struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_def_f5_reduce_to_ext_comparative(int                  j,
                                         int                  k,
                                         int                  l,
                                         struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_f5_reduce_to_ext_comparative(int                  j,
                                     int                  k,
                                     int                  l,
                                     struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_def_f5_reduce_to_stem(int                   j,
                              int                   k,
                              int                   l,
                              struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_user_f5_reduce_to_stem(int                   j,
                          int                   k,
                          int                   l,
                          struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_user_def_f5_reduce_to_stem_comparative(int                   j,
                                          int                   k,
                                          int                   l,
                                          struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_user_f5_reduce_to_stem_comparative(int                   j,
                                      int                   k,
                                      int                   l,
                                      struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_default_f5_split_in_ext_stem(int                   j,
                                int                   k,
                                int                   l,
                                struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_user_def_f5_split_in_ext_stem(int                  j,
                                 int                  k,
                                 int                  l,
                                 struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_f5_split_in_ext_stem(int                  j,
                             int                  k,
                             int                  l,
                             struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_default_f5_split_in_ext_stem_comparative(int                   j,
                                            int                   k,
                                            int                   l,
                                            struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_user_def_f5_split_in_ext_stem_comparative(int                  j,
                                             int                  k,
                                             int                  l,
                                             struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_f5_split_in_ext_stem_comparative(int                  j,
                                         int                  k,
                                         int                  l,
                                         struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_default_f5_split_in_ext_stem1(int                  j,
                                 int                  k,
                                 int                  l,
                                 struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_def_f5_split_in_ext_stem1(int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_user_f5_split_in_ext_stem1(int                   j,
                              int                   k,
                              int                   l,
                              struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_default_f5_split_in_ext_stem1_comparative(int                  j,
                                             int                  k,
                                             int                  l,
                                             struct sc_wrapper_f5 *data);


PRIVATE INLINE int
sc_user_def_f5_split_in_ext_stem1_comparative(int                   j,
                                              int                   k,
                                              int                   l,
                                              struct sc_wrapper_f5  *data);


PRIVATE INLINE int
sc_user_f5_split_in_ext_stem1_comparative(int                   j,
                                          int                   k,
                                          int                   l,
                                          struct sc_wrapper_f5  *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE INLINE void
init_sc_wrapper(vrna_fold_compound_t  *fc,
                struct sc_wrapper_f5  *sc_wrapper)
{
  unsigned int  s;
  vrna_sc_t     *sc, **scs;

  sc_wrapper->up                    = NULL;
  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->n_seq                 = 1;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data_comparative = NULL;

  /* no soft constraints by default */
  sc_wrapper->red_ext       = NULL;
  sc_wrapper->red_stem      = NULL;
  sc_wrapper->decomp_stem   = NULL;
  sc_wrapper->decomp_stem1  = NULL;


  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;
      if (sc) {
        sc_wrapper->up        = sc->energy_up;
        sc_wrapper->user_cb   = sc->f;
        sc_wrapper->user_data = sc->data;

        /* bind correct wrapper functions */
        if (sc->energy_up) {
          if (sc->f) {
            sc_wrapper->red_ext       = &sc_user_def_f5_reduce_to_ext;
            sc_wrapper->red_stem      = &sc_user_def_f5_reduce_to_stem;
            sc_wrapper->decomp_stem   = &sc_user_def_f5_split_in_ext_stem;
            sc_wrapper->decomp_stem1  = &sc_user_def_f5_split_in_ext_stem1;
          } else {
            sc_wrapper->red_ext       = &sc_default_f5_reduce;
            sc_wrapper->red_stem      = &sc_default_f5_reduce;
            sc_wrapper->decomp_stem   = &sc_default_f5_split_in_ext_stem;
            sc_wrapper->decomp_stem1  = &sc_default_f5_split_in_ext_stem1;
          }
        } else if (sc->f) {
          sc_wrapper->red_ext       = &sc_user_f5_reduce_to_ext;
          sc_wrapper->red_stem      = &sc_user_f5_reduce_to_stem;
          sc_wrapper->decomp_stem   = &sc_user_f5_split_in_ext_stem;
          sc_wrapper->decomp_stem1  = &sc_user_f5_split_in_ext_stem1;
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      scs               = fc->scs;
      sc_wrapper->n_seq = fc->n_seq;
      sc_wrapper->a2s   = fc->a2s;
      if (scs) {
        sc_wrapper->up_comparative      = (int ***)vrna_alloc(sizeof(int **) * fc->n_seq);
        sc_wrapper->user_cb_comparative =
          (vrna_callback_sc_energy **)vrna_alloc(sizeof(vrna_callback_sc_energy *) * fc->n_seq);
        sc_wrapper->user_data_comparative = (void **)vrna_alloc(sizeof(void *) * fc->n_seq);

        int provides_sc_up      = 0;
        int provides_sc_user_cb = 0;

        for (s = 0; s < fc->n_seq; s++) {
          if (scs[s]) {
            sc_wrapper->up_comparative[s]         = scs[s]->energy_up;
            sc_wrapper->user_cb_comparative[s]    = scs[s]->f;
            sc_wrapper->user_data_comparative[s]  = scs[s]->data;
            if (scs[s]->energy_up)
              provides_sc_up = 1;

            if (scs[s]->f)
              provides_sc_user_cb = 1;
          }
        }

        /* bind the correct wrapper functions */
        if (provides_sc_up) {
          if (provides_sc_user_cb) {
            sc_wrapper->red_ext       = &sc_user_def_f5_reduce_to_ext_comparative;
            sc_wrapper->red_stem      = &sc_user_def_f5_reduce_to_stem_comparative;
            sc_wrapper->decomp_stem   = &sc_user_def_f5_split_in_ext_stem_comparative;
            sc_wrapper->decomp_stem1  = &sc_user_def_f5_split_in_ext_stem1_comparative;
          } else {
            sc_wrapper->red_ext       = &sc_default_f5_reduce_comparative;
            sc_wrapper->red_stem      = &sc_default_f5_reduce_comparative;
            sc_wrapper->decomp_stem   = &sc_default_f5_split_in_ext_stem_comparative;
            sc_wrapper->decomp_stem1  = &sc_default_f5_split_in_ext_stem1_comparative;
          }
        } else if (provides_sc_user_cb) {
          sc_wrapper->red_ext       = &sc_user_f5_reduce_to_ext_comparative;
          sc_wrapper->red_stem      = &sc_user_f5_reduce_to_stem_comparative;
          sc_wrapper->decomp_stem   = &sc_user_f5_split_in_ext_stem_comparative;
          sc_wrapper->decomp_stem1  = &sc_user_f5_split_in_ext_stem1_comparative;
        }
      }

      break;
  }
}


PRIVATE INLINE void
free_sc_wrapper(struct sc_wrapper_f5 *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}


PRIVATE INLINE int
sc_default_f5_reduce(int                  j,
                     int                  k,
                     int                  l,
                     struct sc_wrapper_f5 *data)
{
  unsigned int  start_2, length_1, length_2;
  int           e_sc, **sc_up;

  sc_up = data->up;

  e_sc = 0;

  length_1  = k - 1;
  start_2   = l + 1;
  length_2  = j - l;

  if (length_1 != 0)
    e_sc += sc_up[1][length_1];

  if (length_2 != 0)
    e_sc += sc_up[start_2][length_2];

  return e_sc;
}


PRIVATE INLINE int
sc_default_f5_reduce_comparative(int                  j,
                                 int                  k,
                                 int                  l,
                                 struct sc_wrapper_f5 *data)
{
  unsigned int  s, start_2, length_1, length_2, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;
  a2s   = data->a2s;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      length_1  = a2s[s][k - 1];
      start_2   = a2s[s][l] + 1;
      length_2  = a2s[s][j] - a2s[s][l];

      if (length_1 != 0)
        e_sc += sc_up[s][1][length_1];

      if (length_2 != 0)
        e_sc += sc_up[s][start_2][length_2];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_user_f5_reduce_to_ext(int                  j,
                         int                  k,
                         int                  l,
                         struct sc_wrapper_f5 *data)
{
  return data->user_cb(1, j, k, l, VRNA_DECOMP_EXT_EXT, data->user_data);
}


PRIVATE INLINE int
sc_user_f5_reduce_to_ext_comparative(int                  j,
                                     int                  k,
                                     int                  l,
                                     struct sc_wrapper_f5 *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc +=
        data->user_cb_comparative[s](1, j, k, l, VRNA_DECOMP_EXT_EXT,
                                     data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_user_def_f5_reduce_to_ext(int                  j,
                             int                  k,
                             int                  l,
                             struct sc_wrapper_f5 *data)
{
  return sc_default_f5_reduce(j, k, l, data) +
         sc_user_f5_reduce_to_ext(j, k, l, data);
}


PRIVATE INLINE int
sc_user_def_f5_reduce_to_ext_comparative(int                  j,
                                         int                  k,
                                         int                  l,
                                         struct sc_wrapper_f5 *data)
{
  return sc_default_f5_reduce_comparative(j, k, l, data) +
         sc_user_f5_reduce_to_ext_comparative(j, k, l, data);
}


PRIVATE INLINE int
sc_user_f5_reduce_to_stem(int                   j,
                          int                   k,
                          int                   l,
                          struct sc_wrapper_f5  *data)
{
  return data->user_cb(1, j, k, l, VRNA_DECOMP_EXT_STEM, data->user_data);
}


PRIVATE INLINE int
sc_user_f5_reduce_to_stem_comparative(int                   j,
                                      int                   k,
                                      int                   l,
                                      struct sc_wrapper_f5  *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc +=
        data->user_cb_comparative[s](1, j, k, l, VRNA_DECOMP_EXT_STEM,
                                     data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_user_def_f5_reduce_to_stem(int                   j,
                              int                   k,
                              int                   l,
                              struct sc_wrapper_f5  *data)
{
  return sc_default_f5_reduce(j, k, l, data) +
         sc_user_f5_reduce_to_stem(j, k, l, data);
}


PRIVATE INLINE int
sc_user_def_f5_reduce_to_stem_comparative(int                   j,
                                          int                   k,
                                          int                   l,
                                          struct sc_wrapper_f5  *data)
{
  return sc_default_f5_reduce_comparative(j, k, l, data) +
         sc_user_f5_reduce_to_stem_comparative(j, k, l, data);
}


PRIVATE INLINE int
sc_default_f5_split_in_ext_stem(int                   j,
                                int                   k,
                                int                   l,
                                struct sc_wrapper_f5  *data)
{
  unsigned int  start_1, length_1;
  int           e_sc, **sc_up;

  sc_up = data->up;

  e_sc = 0;

  start_1   = k + 1;
  length_1  = l - k - 1;

  if (length_1 != 0)
    e_sc += sc_up[start_1][length_1];

  return e_sc;
}


PRIVATE INLINE int
sc_default_f5_split_in_ext_stem_comparative(int                   j,
                                            int                   k,
                                            int                   l,
                                            struct sc_wrapper_f5  *data)
{
  unsigned int  s, start_1, length_1, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;
  a2s   = data->a2s;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][k] + 1;
      length_1  = a2s[l - 1] - a2s[k];

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_user_f5_split_in_ext_stem(int                  j,
                             int                  k,
                             int                  l,
                             struct sc_wrapper_f5 *data)
{
  return data->user_cb(1, j, k, l, VRNA_DECOMP_EXT_EXT_STEM, data->user_data);
}


PRIVATE INLINE int
sc_user_f5_split_in_ext_stem_comparative(int                  j,
                                         int                  k,
                                         int                  l,
                                         struct sc_wrapper_f5 *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc +=
        data->user_cb_comparative[s](1, j, k, l, VRNA_DECOMP_EXT_EXT_STEM,
                                     data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_user_def_f5_split_in_ext_stem(int                  j,
                                 int                  k,
                                 int                  l,
                                 struct sc_wrapper_f5 *data)
{
  return sc_default_f5_split_in_ext_stem(j, k, l, data) +
         sc_user_f5_split_in_ext_stem(j, k, l, data);
}


PRIVATE INLINE int
sc_user_def_f5_split_in_ext_stem_comparative(int                  j,
                                             int                  k,
                                             int                  l,
                                             struct sc_wrapper_f5 *data)
{
  return sc_default_f5_split_in_ext_stem_comparative(j, k, l, data) +
         sc_user_f5_split_in_ext_stem_comparative(j, k, l, data);
}


PRIVATE INLINE int
sc_default_f5_split_in_ext_stem1(int                  j,
                                 int                  k,
                                 int                  l,
                                 struct sc_wrapper_f5 *data)
{
  unsigned int  start_1, start_2, length_1, length_2;
  int           e_sc, **sc_up;

  sc_up = data->up;

  e_sc = 0;

  start_1   = k + 1;
  length_1  = l - k - 1;
  start_2   = j;
  length_2  = 1;

  if (length_1 != 0)
    e_sc += sc_up[start_1][length_1];

  e_sc += sc_up[start_2][length_2];

  return e_sc;
}


PRIVATE INLINE int
sc_default_f5_split_in_ext_stem1_comparative(int                  j,
                                             int                  k,
                                             int                  l,
                                             struct sc_wrapper_f5 *data)
{
  unsigned int  s, start_1, start_2, length_1, length_2, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;

  e_sc = 0;
  a2s  = data->a2s;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][k] + 1;
      length_1  = a2s[s][l - 1] - a2s[s][k];
      start_2   = a2s[s][j];
      length_2  = 1;

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];

      e_sc += sc_up[s][start_2][length_2];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_user_f5_split_in_ext_stem1(int                   j,
                              int                   k,
                              int                   l,
                              struct sc_wrapper_f5  *data)
{
  return data->user_cb(1, j, k, l, VRNA_DECOMP_EXT_EXT_STEM1, data->user_data);
}


PRIVATE INLINE int
sc_user_f5_split_in_ext_stem1_comparative(int                   j,
                                          int                   k,
                                          int                   l,
                                          struct sc_wrapper_f5  *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc +=
        data->user_cb_comparative[s](1, j, k, l, VRNA_DECOMP_EXT_EXT_STEM1,
                                     data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_user_def_f5_split_in_ext_stem1(int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_wrapper_f5  *data)
{
  return sc_default_f5_split_in_ext_stem1(j, k, l, data) +
         sc_user_f5_split_in_ext_stem1(j, k, l, data);
}


PRIVATE INLINE int
sc_user_def_f5_split_in_ext_stem1_comparative(int                   j,
                                              int                   k,
                                              int                   l,
                                              struct sc_wrapper_f5  *data)
{
  return sc_default_f5_split_in_ext_stem1_comparative(j, k, l, data) +
         sc_user_f5_split_in_ext_stem1_comparative(j, k, l, data);
}