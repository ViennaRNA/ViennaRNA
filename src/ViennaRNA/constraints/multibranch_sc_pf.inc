struct sc_mb_exp_dat;

typedef FLT_OR_DBL (*sc_mb_exp_pair_cb)(int                  i,
                                       int                  j,
                                       struct sc_mb_exp_dat *data);


typedef FLT_OR_DBL (*sc_mb_exp_red_cb)(int                   i,
                                      int                   j,
                                      int                   k,
                                      int                   l,
                                      struct sc_mb_exp_dat  *data);


struct sc_mb_exp_dat {
  unsigned int                n;
  unsigned int                n_seq;
  unsigned int                **a2s;

  int                         *idx;

  FLT_OR_DBL                  **up;
  FLT_OR_DBL                  ***up_comparative;
  FLT_OR_DBL                  *bp;
  FLT_OR_DBL                  **bp_comparative;
  FLT_OR_DBL                  **bp_local;
  FLT_OR_DBL                  ***bp_local_comparative;

  sc_mb_exp_pair_cb           pair;
  sc_mb_exp_pair_cb           pair_ext;
  sc_mb_exp_red_cb            red_stem;
  sc_mb_exp_red_cb            red_ml;
  sc_mb_exp_red_cb            decomp_ml;

  vrna_sc_exp_f user_cb;
  void                        *user_data;

  vrna_sc_exp_f *user_cb_comparative;
  void                        **user_data_comparative;
};


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp(int                  i,
                     int                  j,
                     struct sc_mb_exp_dat *data)
{
  return data->bp[data->idx[j] + i];
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_comparative(int                  i,
                                 int                  j,
                                 struct sc_mb_exp_dat *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->bp_comparative[s])
      sc *= data->bp_comparative[s][data->idx[j] + i];

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_local(int                  i,
                           int                  j,
                           struct sc_mb_exp_dat *data)
{
  return data->bp_local[i][j - i];
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_local_comparative(int                  i,
                                       int                  j,
                                       struct sc_mb_exp_dat *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->bp_local_comparative[s])
      sc *= data->bp_local_comparative[s][i][j - i];

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_user(int                  i,
                       int                  j,
                       struct sc_mb_exp_dat *data)
{
  return data->user_cb(i, j, i + 1, j - 1,
                       VRNA_DECOMP_PAIR_ML,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_ext_cb_user(int                  i,
                           int                  j,
                           struct sc_mb_exp_dat *data)
{
  return data->user_cb(i, j, i - 1, j + 1,
                       VRNA_DECOMP_PAIR_ML,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_user_comparative(int                  i,
                                   int                  j,
                                   struct sc_mb_exp_dat *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      sc *= data->user_cb_comparative[s](i, j, i + 1, j - 1,
                                         VRNA_DECOMP_PAIR_ML,
                                         data->user_data_comparative[s]);

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_ext_cb_user_comparative(int                  i,
                                       int                  j,
                                       struct sc_mb_exp_dat *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      sc *= data->user_cb_comparative[s](i, j, i - 1, j + 1,
                                         VRNA_DECOMP_PAIR_ML,
                                         data->user_data_comparative[s]);

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_user(int                   i,
                          int                   j,
                          struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_pair_cb_bp(i, j, data) *
         sc_mb_exp_pair_cb_user(i, j, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_user_comparative(int                   i,
                                      int                   j,
                                      struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_pair_cb_bp_comparative(i, j, data) *
         sc_mb_exp_pair_cb_user_comparative(i, j, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_local_user(int                   i,
                                int                   j,
                                struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_pair_cb_bp_local(i, j, data) *
         sc_mb_exp_pair_cb_user(i, j, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_pair_cb_bp_local_user_comparative(int                   i,
                                            int                   j,
                                            struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_pair_cb_bp_local_comparative(i, j, data) *
         sc_mb_exp_pair_cb_user_comparative(i, j, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_up(int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    struct sc_mb_exp_dat  *data)
{
  int         l1  = k - i;
  int         l2  = j - l;
  FLT_OR_DBL  sc  = 1.;

  if (l1 > 0)
    sc *= data->up[i][l1];

  if (l2 > 0)
    sc *= data->up[l + 1][l2];

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_up_comparative(int                   i,
                                int                   j,
                                int                   k,
                                int                   l,
                                struct sc_mb_exp_dat  *data)
{
  unsigned int  s, l1, l2, start1, start2;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++) {
    if (data->up_comparative[s]) {
      start1  = data->a2s[s][i];
      l1      = data->a2s[s][k] - data->a2s[s][i];
      start2  = data->a2s[s][l] + 1;
      l2      = data->a2s[s][j] - data->a2s[s][l];

      if (l1 != 0)
        sc *= data->up_comparative[s][start1][l1];

      if (l2 != 0)
        sc *= data->up_comparative[s][start2][l2];
    }
  }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_user(int                   i,
                      int                   j,
                      int                   k,
                      int                   l,
                      struct sc_mb_exp_dat  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_ML,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_user_comparative(int                   i,
                                  int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_mb_exp_dat  *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      sc *= data->user_cb_comparative[s](i, j, k, l,
                                         VRNA_DECOMP_ML_ML,
                                         data->user_data_comparative[s]);

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_up_user(int                  i,
                         int                  j,
                         int                  k,
                         int                  l,
                         struct sc_mb_exp_dat *data)
{
  return sc_mb_exp_red_cb_up(i, j, k, l, data) *
         sc_mb_exp_red_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_up_user_comparative(int                  i,
                                     int                  j,
                                     int                  k,
                                     int                  l,
                                     struct sc_mb_exp_dat *data)
{
  return sc_mb_exp_red_cb_up_comparative(i, j, k, l, data) *
         sc_mb_exp_red_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_stem_user(int                  i,
                           int                  j,
                           int                  k,
                           int                  l,
                           struct sc_mb_exp_dat *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_STEM,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_stem_user_comparative(int                  i,
                                       int                  j,
                                       int                  k,
                                       int                  l,
                                       struct sc_mb_exp_dat *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      sc *= data->user_cb_comparative[s](i, j, k, l,
                                         VRNA_DECOMP_ML_STEM,
                                         data->user_data);

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_stem_up_user(int                   i,
                              int                   j,
                              int                   k,
                              int                   l,
                              struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_red_cb_up(i, j, k, l, data) *
         sc_mb_exp_red_cb_stem_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_red_cb_stem_up_user_comparative(int                   i,
                                          int                   j,
                                          int                   k,
                                          int                   l,
                                          struct sc_mb_exp_dat  *data)
{
  return sc_mb_exp_red_cb_up_comparative(i, j, k, l, data) *
         sc_mb_exp_red_cb_stem_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_split_cb_user(int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        struct sc_mb_exp_dat  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_ML_ML,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_mb_exp_split_cb_user_comparative(int                   i,
                                    int                   j,
                                    int                   k,
                                    int                   l,
                                    struct sc_mb_exp_dat  *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      sc *= data->user_cb_comparative[s](i, j, k, l,
                                         VRNA_DECOMP_ML_ML_ML,
                                         data->user_data_comparative[s]);

  return sc;
}


PRIVATE INLINE void
init_sc_mb_exp(vrna_fold_compound_t *fc,
               struct sc_mb_exp_dat *sc_wrapper)
{
  unsigned char sliding_window;
  vrna_sc_t     *sc, **scs;

  sc_wrapper->n     = fc->length;
  sc_wrapper->n_seq = 1;
  sc_wrapper->idx   = fc->jindx;
  sc_wrapper->a2s   = NULL;

  sc_wrapper->up                    = NULL;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->bp                    = NULL;
  sc_wrapper->bp_comparative        = NULL;
  sc_wrapper->bp_local              = NULL;
  sc_wrapper->bp_local_comparative  = NULL;

  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data_comparative = NULL;

  sc_wrapper->pair      = NULL;
  sc_wrapper->pair_ext  = NULL;
  sc_wrapper->red_stem  = NULL;
  sc_wrapper->red_ml    = NULL;
  sc_wrapper->decomp_ml = NULL;

  sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        unsigned int provides_sc_up, provides_sc_bp, provides_sc_user;

        provides_sc_up    = 0;
        provides_sc_bp    = 0;
        provides_sc_user  = 0;

        sc_wrapper->up        = sc->exp_energy_up;
        sc_wrapper->user_cb   = sc->exp_f;
        sc_wrapper->user_data = sc->data;

        if (sliding_window)
          sc_wrapper->bp_local = sc->exp_energy_bp_local;
        else
          sc_wrapper->bp = sc->exp_energy_bp;

        if (sc->exp_energy_up)
          provides_sc_up = 1;

        if (sliding_window) {
          if (sc->exp_energy_bp_local)
            provides_sc_bp = 1;
        } else if (sc->exp_energy_bp) {
          provides_sc_bp = 1;
        }

        if (sc->exp_f)
          provides_sc_user = 1;

        /* done initializing, now assign function pointers */
        if (provides_sc_user) {
          sc_wrapper->decomp_ml = &sc_mb_exp_split_cb_user;
          sc_wrapper->red_stem  = &sc_mb_exp_red_cb_stem_user;
          sc_wrapper->red_ml    = &sc_mb_exp_red_cb_user;
          sc_wrapper->pair      = &sc_mb_exp_pair_cb_user;
          if (!sliding_window)
            sc_wrapper->pair_ext  = &sc_mb_exp_pair_ext_cb_user;

          if (provides_sc_bp) {
            if (sliding_window) {
              sc_wrapper->pair = &sc_mb_exp_pair_cb_bp_local_user;
            } else {
              sc_wrapper->pair      = &sc_mb_exp_pair_cb_bp_user;
              sc_wrapper->pair_ext  = &sc_mb_exp_pair_ext_cb_user;
            }
          }

          if (provides_sc_up) {
            sc_wrapper->red_stem  = &sc_mb_exp_red_cb_stem_up_user;
            sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up_user;
          }
        } else if (provides_sc_bp) {
          if (sliding_window) {
            sc_wrapper->pair = &sc_mb_exp_pair_cb_bp_local;
          } else {
            sc_wrapper->pair      = &sc_mb_exp_pair_cb_bp;
          }

          if (provides_sc_up) {
            sc_wrapper->red_stem  = &sc_mb_exp_red_cb_up;
            sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up;
          }
        } else if (provides_sc_up) {
          sc_wrapper->red_stem  = &sc_mb_exp_red_cb_up;
          sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up;
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      sc_wrapper->a2s   = fc->a2s;
      sc_wrapper->n_seq = fc->n_seq;
      scs               = fc->scs;

      if (scs) {
        unsigned int s, provides_sc_up, provides_sc_bp, provides_sc_user;

        provides_sc_up    = 0;
        provides_sc_bp    = 0;
        provides_sc_user  = 0;

        sc_wrapper->up_comparative = (FLT_OR_DBL ***)vrna_alloc(
          sizeof(FLT_OR_DBL * *) *
          fc->n_seq);

        sc_wrapper->bp_comparative = (FLT_OR_DBL **)vrna_alloc(
          sizeof(FLT_OR_DBL *) *
          fc->n_seq);

        sc_wrapper->bp_local_comparative = (FLT_OR_DBL ***)vrna_alloc(
          sizeof(FLT_OR_DBL * *) *
          fc->n_seq);

        sc_wrapper->user_cb_comparative = (vrna_sc_exp_f *)vrna_alloc(
          sizeof(vrna_sc_exp_f) *
          fc->n_seq);

        sc_wrapper->user_data_comparative = (void **)vrna_alloc(
          sizeof(void *) *
          fc->n_seq);

        for (s = 0; s < fc->n_seq; s++) {
          if (scs[s]) {
            sc_wrapper->up_comparative[s]       = scs[s]->exp_energy_up;
            sc_wrapper->bp_comparative[s]       = (sliding_window) ? NULL : scs[s]->exp_energy_bp;
            sc_wrapper->bp_local_comparative[s] =
              (sliding_window) ? scs[s]->exp_energy_bp_local : NULL;
            sc_wrapper->user_cb_comparative[s]    = scs[s]->exp_f;
            sc_wrapper->user_data_comparative[s]  = scs[s]->data;

            if (scs[s]->exp_energy_up)
              provides_sc_up = 1;

            if (sliding_window) {
              if (scs[s]->exp_energy_bp_local)
                provides_sc_bp = 1;
            } else if (scs[s]->exp_energy_bp) {
              provides_sc_bp = 1;
            }

            if (scs[s]->exp_f)
              provides_sc_user = 1;
          }
        }

        /* done initializing, now assign function pointers */
        if (provides_sc_user) {
          sc_wrapper->decomp_ml = &sc_mb_exp_split_cb_user_comparative;
          sc_wrapper->red_stem  = &sc_mb_exp_red_cb_stem_user_comparative;
          sc_wrapper->red_ml    = &sc_mb_exp_red_cb_user_comparative;
          sc_wrapper->pair      = &sc_mb_exp_pair_cb_user_comparative;
          if (!sliding_window)
            sc_wrapper->pair_ext  = &sc_mb_exp_pair_ext_cb_user_comparative;

          if (provides_sc_bp) {
            if (sliding_window) {
              sc_wrapper->pair = &sc_mb_exp_pair_cb_bp_local_user_comparative;
            } else {
              sc_wrapper->pair      = &sc_mb_exp_pair_cb_bp_user_comparative;
              sc_wrapper->pair_ext  = &sc_mb_exp_pair_ext_cb_user_comparative;
            }
          }

          if (provides_sc_up) {
            sc_wrapper->red_stem  = &sc_mb_exp_red_cb_stem_up_user_comparative;
            sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up_user_comparative;
          }
        } else if (provides_sc_bp) {
          if (sliding_window) {
            sc_wrapper->pair = &sc_mb_exp_pair_cb_bp_local_comparative;
          } else {
            sc_wrapper->pair      = &sc_mb_exp_pair_cb_bp_comparative;
          }

          if (provides_sc_up) {
            sc_wrapper->red_stem  = &sc_mb_exp_red_cb_up_comparative;
            sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up_comparative;
          }
        } else if (provides_sc_up) {
          sc_wrapper->red_stem  = &sc_mb_exp_red_cb_up_comparative;
          sc_wrapper->red_ml    = &sc_mb_exp_red_cb_up_comparative;
        }
      }

      break;
  }
}


PRIVATE INLINE void
free_sc_mb_exp(struct sc_mb_exp_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->bp_comparative);
  free(sc_wrapper->bp_local_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}
