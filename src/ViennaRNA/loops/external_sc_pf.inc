/*
 *  This file contains soft constraints wrappers for exterior loop evaluation (partition function version)
 */
struct sc_wrapper_exp_ext;

typedef FLT_OR_DBL (sc_ext_exp_red)(int                       i,
                                    int                       j,
                                    int                       k,
                                    int                       l,
                                    struct sc_wrapper_exp_ext *data);
typedef FLT_OR_DBL (sc_ext_exp_stem)(int                        i,
                                     int                        j,
                                     int                        k,
                                     int                        l,
                                     struct sc_wrapper_exp_ext  *data);
typedef FLT_OR_DBL (sc_ext_exp_red_up)(int                        i,
                                       int                        j,
                                       struct sc_wrapper_exp_ext  *data);
typedef FLT_OR_DBL (sc_ext_exp_split)(int                       i,
                                      int                       j,
                                      int                       k,
                                      struct sc_wrapper_exp_ext *data);

struct sc_wrapper_exp_ext {
  FLT_OR_DBL                  **up;

  sc_ext_exp_red              *red_ext;
  sc_ext_exp_stem             *red_stem;
  sc_ext_exp_red_up           *red_up;
  sc_ext_exp_split            *split;

  vrna_callback_sc_exp_energy *user_cb;
  void                        *user_data;

  /* below attributes are for comparative structure prediction */
  int                         n_seq;
  unsigned int                **a2s;
  FLT_OR_DBL                  ***up_comparative;

  vrna_callback_sc_exp_energy **user_cb_comparative;
  void                        **user_data_comparative;
};


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_reduce(int                       i,
                  int                       j,
                  int                       k,
                  int                       l,
                  struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_reduce_comparative(int                       i,
                              int                       j,
                              int                       k,
                              int                       l,
                              struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_reduce_to_ext(int                       i,
                              int                       j,
                              int                       k,
                              int                       l,
                              struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_reduce_to_ext_comparative(int                       i,
                                          int                       j,
                                          int                       k,
                                          int                       l,
                                          struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_reduce_to_ext(int                       i,
                                  int                       j,
                                  int                       k,
                                  int                       l,
                                  struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_reduce_to_ext_comparative(int                       i,
                                              int                       j,
                                              int                       k,
                                              int                       l,
                                              struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_reduce_to_stem(int                        i,
                               int                        j,
                               int                        k,
                               int                        l,
                               struct sc_wrapper_exp_ext  *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_reduce_to_stem_comparative(int                        i,
                                           int                        j,
                                           int                        k,
                                           int                        l,
                                           struct sc_wrapper_exp_ext  *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_reduce_to_stem(int                        i,
                                   int                        j,
                                   int                        k,
                                   int                        l,
                                   struct sc_wrapper_exp_ext  *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_reduce_to_stem_comparative(int                        i,
                                               int                        j,
                                               int                        k,
                                               int                        l,
                                               struct sc_wrapper_exp_ext  *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_up(int                       i,
              int                       j,
              struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_up_comparative(int                       i,
                          int                       j,
                          struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_up(int                        i,
                   int                        j,
                   struct sc_wrapper_exp_ext  *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_up_comparative(int                        i,
                               int                        j,
                               struct sc_wrapper_exp_ext  *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_up(int                        i,
                       int                        j,
                       struct sc_wrapper_exp_ext  *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_up_comparative(int                        i,
                                   int                        j,
                                   struct sc_wrapper_exp_ext  *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_split(int                       i,
                      int                       j,
                      int                       k,
                      struct sc_wrapper_exp_ext *data);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_split_comparative(int                       i,
                                  int                       j,
                                  int                       k,
                                  struct sc_wrapper_exp_ext *data);


PRIVATE INLINE void
init_sc_wrapper_ext(vrna_fold_compound_t      *fc,
                    struct sc_wrapper_exp_ext *sc_wrapper);


PRIVATE INLINE void
free_sc_wrapper_ext(struct sc_wrapper_exp_ext *sc_wrapper);


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_reduce(int                       i,
                  int                       j,
                  int                       k,
                  int                       l,
                  struct sc_wrapper_exp_ext *data)
{
  unsigned int  start_2, length_1, length_2;
  FLT_OR_DBL    q_sc, **sc_up;

  sc_up = data->up;

  q_sc = 1.;

  length_1  = k - i;
  start_2   = l + 1;
  length_2  = j - l;

  if (length_1 != 0)
    q_sc *= sc_up[i][length_1];

  if (length_2 != 0)
    q_sc *= sc_up[start_2][length_2];

  return q_sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_reduce_comparative(int                       i,
                              int                       j,
                              int                       k,
                              int                       l,
                              struct sc_wrapper_exp_ext *data)
{
  unsigned int  s, start_2, length_1, length_2, **a2s;
  FLT_OR_DBL    q_sc, ***sc_up;

  sc_up = data->up_comparative;
  a2s   = data->a2s;

  q_sc = 1.;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      length_1  = a2s[s][k] - a2s[s][i];
      start_2   = a2s[s][l] + 1;
      length_2  = a2s[s][j] - a2s[s][l];

      if (length_1 != 0)
        q_sc *= sc_up[s][a2s[s][i]][length_1];

      if (length_2 != 0)
        q_sc *= sc_up[s][start_2][length_2];
    }
  }

  return q_sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_reduce_to_ext(int                       i,
                              int                       j,
                              int                       k,
                              int                       l,
                              struct sc_wrapper_exp_ext *data)
{
  return data->user_cb(i, j, k, l, VRNA_DECOMP_EXT_EXT, data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_reduce_to_ext_comparative(int                       i,
                                          int                       j,
                                          int                       k,
                                          int                       l,
                                          struct sc_wrapper_exp_ext *data)
{
  unsigned int  s;
  FLT_OR_DBL    q_sc;

  q_sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    q_sc *=
      data->user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_EXT_EXT, data->user_data_comparative[s]);

  return q_sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_reduce_to_ext(int                       i,
                                  int                       j,
                                  int                       k,
                                  int                       l,
                                  struct sc_wrapper_exp_ext *data)
{
  return sc_ext_exp_reduce(i, j, k, l, data) +
         sc_ext_exp_user_reduce_to_ext(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_reduce_to_ext_comparative(int                       i,
                                              int                       j,
                                              int                       k,
                                              int                       l,
                                              struct sc_wrapper_exp_ext *data)
{
  return sc_ext_exp_reduce_comparative(i, j, k, l, data) +
         sc_ext_exp_user_reduce_to_ext_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_reduce_to_stem(int                        i,
                               int                        j,
                               int                        k,
                               int                        l,
                               struct sc_wrapper_exp_ext  *data)
{
  return data->user_cb(i, j, k, l, VRNA_DECOMP_EXT_STEM, data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_reduce_to_stem_comparative(int                        i,
                                           int                        j,
                                           int                        k,
                                           int                        l,
                                           struct sc_wrapper_exp_ext  *data)
{
  unsigned int  s;
  FLT_OR_DBL    q_sc;

  q_sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    q_sc *=
      data->user_cb_comparative[s](i, j, k, l, VRNA_DECOMP_EXT_STEM,
                                   data->user_data_comparative[s]);

  return q_sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_reduce_to_stem(int                        i,
                                   int                        j,
                                   int                        k,
                                   int                        l,
                                   struct sc_wrapper_exp_ext  *data)
{
  return sc_ext_exp_reduce(i, j, k, l, data) +
         sc_ext_exp_user_reduce_to_stem(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_reduce_to_stem_comparative(int                        i,
                                               int                        j,
                                               int                        k,
                                               int                        l,
                                               struct sc_wrapper_exp_ext  *data)
{
  return sc_ext_exp_reduce_comparative(i, j, k, l, data) +
         sc_ext_exp_user_reduce_to_stem_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_up(int                       i,
              int                       j,
              struct sc_wrapper_exp_ext *data)
{
  unsigned int  length;
  FLT_OR_DBL    q_sc, **sc_up;

  sc_up   = data->up;
  length  = j - i + 1;
  q_sc    = 1.;

  if (length != 0)
    q_sc *= sc_up[i][length];

  return q_sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_up_comparative(int                       i,
                          int                       j,
                          struct sc_wrapper_exp_ext *data)
{
  unsigned int  length, s, **a2s;
  FLT_OR_DBL    q_sc, ***sc_up;

  a2s   = data->a2s;
  sc_up = data->up_comparative;

  q_sc = 1.;

  for (s = 0; s < data->n_seq; s++) {
    length = a2s[s][j - 1] - a2s[s][i];

    if (length != 0)
      q_sc *= sc_up[s][a2s[s][i]][length];
  }

  return q_sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_up(int                        i,
                   int                        j,
                   struct sc_wrapper_exp_ext  *data)
{
  return data->user_cb(i, j, i, j, VRNA_DECOMP_EXT_UP, data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_up_comparative(int                        i,
                               int                        j,
                               struct sc_wrapper_exp_ext  *data)
{
  unsigned int  s;
  FLT_OR_DBL    q_sc;

  q_sc = 1.;
  for (s = 0; s < data->n_seq; s++)
    q_sc *=
      data->user_cb_comparative[s](i, j, i, j, VRNA_DECOMP_EXT_UP, data->user_data_comparative[s]);

  return q_sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_up(int                        i,
                       int                        j,
                       struct sc_wrapper_exp_ext  *data)
{
  return sc_ext_exp_up(i, j, data) +
         sc_ext_exp_user_up(i, j, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_def_up_comparative(int                        i,
                                   int                        j,
                                   struct sc_wrapper_exp_ext  *data)
{
  return sc_ext_exp_up_comparative(i, j, data) +
         sc_ext_exp_user_up_comparative(i, j, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_split(int                       i,
                      int                       j,
                      int                       k,
                      struct sc_wrapper_exp_ext *data)
{
  return data->user_cb(i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT, data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_ext_exp_user_split_comparative(int                       i,
                                  int                       j,
                                  int                       k,
                                  struct sc_wrapper_exp_ext *data)
{
  unsigned int  s;
  FLT_OR_DBL    q_sc;

  q_sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    q_sc *=
      data->user_cb_comparative[s](i, j, k - 1, k, VRNA_DECOMP_EXT_EXT_EXT,
                                   data->user_data_comparative[s]);

  return q_sc;
}


PRIVATE INLINE void
init_sc_wrapper_ext(vrna_fold_compound_t      *fc,
                    struct sc_wrapper_exp_ext *sc_wrapper)
{
  vrna_sc_t *sc, **scs;

  sc_wrapper->up                    = NULL;
  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->n_seq                 = 1;
  sc_wrapper->a2s                   = NULL;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data_comparative = NULL;

  /* no soft constraints by default */
  sc_wrapper->red_ext   = NULL;
  sc_wrapper->red_stem  = NULL;
  sc_wrapper->red_up    = NULL;
  sc_wrapper->split     = NULL;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        sc_wrapper->up        = sc->exp_energy_up;
        sc_wrapper->user_cb   = sc->exp_f;
        sc_wrapper->user_data = sc->data;

        /* bind correct wrapper functions */
        if (sc->exp_energy_up) {
          if (sc->exp_f) {
            sc_wrapper->red_ext   = &sc_ext_exp_user_def_reduce_to_ext;
            sc_wrapper->red_stem  = &sc_ext_exp_user_def_reduce_to_stem;
            sc_wrapper->red_up    = &sc_ext_exp_user_def_up;
            sc_wrapper->split     = &sc_ext_exp_user_split;
          } else {
            sc_wrapper->red_ext   = &sc_ext_exp_reduce;
            sc_wrapper->red_stem  = &sc_ext_exp_reduce;
            sc_wrapper->red_up    = &sc_ext_exp_up;
          }
        } else if (sc->exp_f) {
          sc_wrapper->red_ext   = &sc_ext_exp_user_reduce_to_ext;
          sc_wrapper->red_stem  = &sc_ext_exp_user_reduce_to_stem;
          sc_wrapper->red_up    = &sc_ext_exp_user_up;
          sc_wrapper->split     = &sc_ext_exp_user_split;
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      scs               = fc->scs;
      sc_wrapper->n_seq = fc->n_seq;
      sc_wrapper->a2s   = fc->a2s;

      if (scs) {
        unsigned int s;

        sc_wrapper->up_comparative = (FLT_OR_DBL ***)vrna_alloc(
          sizeof(FLT_OR_DBL * *) * fc->n_seq);
        sc_wrapper->user_cb_comparative =
          (vrna_callback_sc_exp_energy **)vrna_alloc(
            sizeof(vrna_callback_sc_exp_energy *) * fc->n_seq);
        sc_wrapper->user_data_comparative = (void **)vrna_alloc(sizeof(void *) * fc->n_seq);

        int provides_sc_up      = 0;
        int provides_sc_user_cb = 0;

        for (s = 0; s < fc->n_seq; s++) {
          if (scs[s]) {
            sc_wrapper->up_comparative[s]         = scs[s]->exp_energy_up;
            sc_wrapper->user_cb_comparative[s]    = scs[s]->exp_f;
            sc_wrapper->user_data_comparative[s]  = scs[s]->data;
            if (scs[s]->exp_energy_up)
              provides_sc_up = 1;

            if (scs[s]->exp_f)
              provides_sc_user_cb = 1;

            /* bind the correct wrapper functions */
            if (provides_sc_up) {
              if (provides_sc_user_cb) {
                sc_wrapper->red_ext   = &sc_ext_exp_user_def_reduce_to_ext_comparative;
                sc_wrapper->red_stem  = &sc_ext_exp_user_def_reduce_to_stem_comparative;
                sc_wrapper->red_up    = &sc_ext_exp_user_def_up_comparative;
                sc_wrapper->split     = &sc_ext_exp_user_split_comparative;
              } else {
                sc_wrapper->red_ext   = &sc_ext_exp_reduce_comparative;
                sc_wrapper->red_stem  = &sc_ext_exp_reduce_comparative;
                sc_wrapper->red_up    = &sc_ext_exp_up_comparative;
              }
            } else if (provides_sc_user_cb) {
              sc_wrapper->red_ext   = &sc_ext_exp_user_reduce_to_ext_comparative;
              sc_wrapper->red_stem  = &sc_ext_exp_user_reduce_to_stem_comparative;
              sc_wrapper->red_up    = &sc_ext_exp_user_up_comparative;
              sc_wrapper->split     = &sc_ext_exp_user_split_comparative;
            }
          }
        }
      }

      break;
  }
}


PRIVATE INLINE void
free_sc_wrapper_ext(struct sc_wrapper_exp_ext *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}
