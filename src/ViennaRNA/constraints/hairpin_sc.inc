struct sc_hp_dat;

typedef int (*sc_hp_cb)(int              i,
                       int              j,
                       struct sc_hp_dat *data);


struct sc_hp_dat {
  unsigned int            n;
  unsigned int            n_seq;
  unsigned int            **a2s;
  int                     *idx;

  int                     **up;
  int                     ***up_comparative;
  int                     *bp;
  int                     **bp_comparative;
  int                     **bp_local;
  int                     ***bp_local_comparative;

  vrna_sc_f user_cb;
  void                    *user_data;

  vrna_sc_f *user_cb_comparative;
  void                    **user_data_comparative;

  sc_hp_cb                pair;
  sc_hp_cb                pair_ext;
};


PRIVATE INLINE int
sc_hp_cb_up(int               i,
            int               j,
            struct sc_hp_dat  *data)
{
  return data->up[i + 1][j - i - 1];
}


PRIVATE INLINE int
sc_hp_cb_up_comparative(int               i,
                        int               j,
                        struct sc_hp_dat  *data)
{
  unsigned int  s, u;
  int           sc;

  sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->up_comparative[s]) {
      u   = data->a2s[s][j - 1] - data->a2s[s][i];
      if (u != 0)
        sc += data->up_comparative[s][data->a2s[s][i] + 1][u];
    }

  return sc;
}


PRIVATE INLINE int
sc_hp_cb_bp(int               i,
            int               j,
            struct sc_hp_dat  *data)
{
  return data->bp[data->idx[j] + i];
}


PRIVATE INLINE int
sc_hp_cb_bp_comparative(int               i,
                        int               j,
                        struct sc_hp_dat  *data)
{
  unsigned int  s;
  int           sc;

  sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->bp_comparative[s])
      sc += data->bp_comparative[s][data->idx[j] + i];

  return sc;
}


PRIVATE INLINE int
sc_hp_cb_bp_local(int               i,
                  int               j,
                  struct sc_hp_dat  *data)
{
  return data->bp_local[i][j - i];
}


PRIVATE INLINE int
sc_hp_cb_bp_local_comparative(int               i,
                              int               j,
                              struct sc_hp_dat  *data)
{
  unsigned int  s;
  int           sc;

  sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->bp_local_comparative[s])
      sc += data->bp_local_comparative[s][i][j - i];

  return sc;
}


PRIVATE INLINE int
sc_hp_cb_user(int               i,
              int               j,
              struct sc_hp_dat  *data)
{
  return data->user_cb(i, j, i, j,
                       VRNA_DECOMP_PAIR_HP,
                       data->user_data);
}


PRIVATE INLINE int
sc_hp_cb_user_comparative(int               i,
                          int               j,
                          struct sc_hp_dat  *data)
{
  unsigned int  s;
  int           sc;

  sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      sc += data->user_cb_comparative[s](i, j, i, j,
                                         VRNA_DECOMP_PAIR_HP,
                                         data->user_data_comparative[s]);

  return sc;
}


PRIVATE INLINE int
sc_hp_cb_up_bp(int              i,
               int              j,
               struct sc_hp_dat *data)
{
  return sc_hp_cb_up(i, j, data) +
         sc_hp_cb_bp(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_bp_comparative(int              i,
                           int              j,
                           struct sc_hp_dat *data)
{
  return sc_hp_cb_up_comparative(i, j, data) +
         sc_hp_cb_bp_comparative(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_bp_local(int              i,
                     int              j,
                     struct sc_hp_dat *data)
{
  return sc_hp_cb_up(i, j, data) +
         sc_hp_cb_bp_local(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_bp_local_comparative(int              i,
                                 int              j,
                                 struct sc_hp_dat *data)
{
  return sc_hp_cb_up_comparative(i, j, data) +
         sc_hp_cb_bp_local_comparative(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_user(int              i,
                 int              j,
                 struct sc_hp_dat *data)
{
  return sc_hp_cb_up(i, j, data) +
         sc_hp_cb_user(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_user_comparative(int              i,
                             int              j,
                             struct sc_hp_dat *data)
{
  return sc_hp_cb_up_comparative(i, j, data) +
         sc_hp_cb_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_bp_user(int              i,
                 int              j,
                 struct sc_hp_dat *data)
{
  return sc_hp_cb_bp(i, j, data) +
         sc_hp_cb_user(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_bp_user_comparative(int              i,
                             int              j,
                             struct sc_hp_dat *data)
{
  return sc_hp_cb_bp_comparative(i, j, data) +
         sc_hp_cb_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_bp_local_user(int              i,
                       int              j,
                       struct sc_hp_dat *data)
{
  return sc_hp_cb_bp_local(i, j, data) +
         sc_hp_cb_user(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_bp_local_user_comparative(int              i,
                                   int              j,
                                   struct sc_hp_dat *data)
{
  return sc_hp_cb_bp_local_comparative(i, j, data) +
         sc_hp_cb_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_bp_user(int               i,
                    int               j,
                    struct sc_hp_dat  *data)
{
  return sc_hp_cb_up(i, j, data) +
         sc_hp_cb_bp(i, j, data) +
         sc_hp_cb_user(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_bp_user_comparative(int               i,
                                int               j,
                                struct sc_hp_dat  *data)
{
  return sc_hp_cb_up_comparative(i, j, data) +
         sc_hp_cb_bp_comparative(i, j, data) +
         sc_hp_cb_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_bp_local_user(int               i,
                          int               j,
                          struct sc_hp_dat  *data)
{
  return sc_hp_cb_up(i, j, data) +
         sc_hp_cb_bp_local(i, j, data) +
         sc_hp_cb_user(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_up_bp_local_user_comparative(int               i,
                                      int               j,
                                      struct sc_hp_dat  *data)
{
  return sc_hp_cb_up_comparative(i, j, data) +
         sc_hp_cb_bp_local_comparative(i, j, data) +
         sc_hp_cb_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_ext_up(int               i,
                int               j,
                struct sc_hp_dat  *data)
{
  int sc, u1, u2;

  u1  = data->n - j;
  u2  = i - 1;
  sc  = 0;

  if (u1 > 0)
    sc += data->up[j + 1][u1];

  if (u2 > 0)
    sc += data->up[1][u2];

  return sc;
}


PRIVATE INLINE int
sc_hp_cb_ext_up_comparative(int               i,
                            int               j,
                            struct sc_hp_dat  *data)
{
  unsigned int  s, u1, u2;
  int           sc;

  sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->up_comparative[s]) {
      u1  = data->a2s[s][data->n] - data->a2s[s][j];
      u2  = (i > 1) ? data->a2s[s][i - 1] - data->a2s[s][1] : 0;

      if (u1 != 0)
        sc += data->up[data->a2s[s][j] + 1][u1];

      if (u2 != 0)
        sc += data->up[data->a2s[s][1]][u2];
    }

  return sc;
}


PRIVATE INLINE int
sc_hp_cb_ext_user(int               i,
                  int               j,
                  struct sc_hp_dat  *data)
{
  return data->user_cb(j, i, j, i,
                       VRNA_DECOMP_PAIR_HP,
                       data->user_data);
}


PRIVATE INLINE int
sc_hp_cb_ext_user_comparative(int               i,
                              int               j,
                              struct sc_hp_dat  *data)
{
  unsigned int  s;
  int           sc;

  sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      sc += data->user_cb_comparative[s](j, i, j, i,
                                         VRNA_DECOMP_PAIR_HP,
                                         data->user_data_comparative[s]);

  return sc;
}


PRIVATE INLINE int
sc_hp_cb_ext_up_user(int              i,
                     int              j,
                     struct sc_hp_dat *data)
{
  return sc_hp_cb_ext_up(i, j, data) +
         sc_hp_cb_ext_user(i, j, data);
}


PRIVATE INLINE int
sc_hp_cb_ext_up_user_comparative(int              i,
                                 int              j,
                                 struct sc_hp_dat *data)
{
  return sc_hp_cb_ext_up_comparative(i, j, data) +
         sc_hp_cb_ext_user_comparative(i, j, data);
}


PRIVATE INLINE void
init_sc_hp(vrna_fold_compound_t *fc,
           struct sc_hp_dat     *sc_wrapper)
{
  unsigned char sliding_window;
  unsigned int  s;
  vrna_sc_t     *sc, **scs;

  if (fc->matrices)
    sliding_window = (fc->matrices->type == VRNA_MX_WINDOW) ? 1 : 0;
  else if ((fc->type == VRNA_FC_TYPE_SINGLE) && (fc->sc))
    sliding_window = (fc->sc->type == VRNA_SC_WINDOW) ? 1 : 0;
  else if (fc->hc)
    sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  else
    sliding_window = 0;

  sc_wrapper->n     = (int)fc->length;
  sc_wrapper->idx   = fc->jindx;
  sc_wrapper->n_seq = 1;
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

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        unsigned int provides_sc_up, provides_sc_bp, provides_sc_user;

        provides_sc_up    = 0;
        provides_sc_bp    = 0;
        provides_sc_user  = 0;

        sc_wrapper->up        = sc->energy_up;
        sc_wrapper->bp        = (sliding_window) ? NULL : sc->energy_bp;
        sc_wrapper->bp_local  = (sliding_window) ? sc->energy_bp_local : NULL;
        sc_wrapper->user_cb   = sc->f;
        sc_wrapper->user_data = sc->data;

        if (sc->energy_up)
          provides_sc_up = 1;

        if (sliding_window) {
          if (sc->energy_bp_local)
            provides_sc_bp = 1;
        } else if (sc->energy_bp) {
          provides_sc_bp = 1;
        }

        if (sc->f)
          provides_sc_user = 1;

        if (provides_sc_user) {
          sc_wrapper->pair_ext = &sc_hp_cb_ext_user;
          if (provides_sc_up) {
            sc_wrapper->pair_ext = &sc_hp_cb_ext_up_user;

            if (provides_sc_bp)
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_hp_cb_up_bp_local_user :
                                 &sc_hp_cb_up_bp_user;
            else
              sc_wrapper->pair = &sc_hp_cb_up_user;
          } else if (provides_sc_bp) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_hp_cb_bp_local_user :
                               &sc_hp_cb_bp_user;
          } else {
            sc_wrapper->pair = &sc_hp_cb_user;
          }
        } else {
          if (provides_sc_up) {
            sc_wrapper->pair_ext = &sc_hp_cb_ext_up;
            if (provides_sc_bp) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_hp_cb_up_bp_local :
                                 &sc_hp_cb_up_bp;
            } else {
              sc_wrapper->pair = &sc_hp_cb_up;
            }
          } else if (provides_sc_bp) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_hp_cb_bp_local :
                               &sc_hp_cb_bp;
          }
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      sc_wrapper->n_seq = fc->n_seq;
      sc_wrapper->a2s   = fc->a2s;

      scs = fc->scs;

      if (scs) {
        unsigned int provides_sc_up, provides_sc_bp, provides_sc_user;

        provides_sc_up    = 0;
        provides_sc_bp    = 0;
        provides_sc_user  = 0;

        sc_wrapper->up_comparative = (int ***)vrna_alloc(
          sizeof(int **) *
          fc->n_seq);

        sc_wrapper->bp_comparative = (int **)vrna_alloc(
          sizeof(int *) *
          fc->n_seq);

        sc_wrapper->bp_local_comparative = (int ***)vrna_alloc(
          sizeof(int **) *
          fc->n_seq);

        sc_wrapper->user_cb_comparative = (vrna_sc_f *)vrna_alloc(
          sizeof(vrna_sc_f) *
          fc->n_seq);

        sc_wrapper->user_data_comparative = (void **)vrna_alloc(
          sizeof(void *) *
          fc->n_seq);

        for (s = 0; s < fc->n_seq; s++) {
          if (scs[s]) {
            sliding_window                      = (scs[s]->type == VRNA_SC_WINDOW) ? 1 : 0;
            sc_wrapper->up_comparative[s]       = scs[s]->energy_up;
            sc_wrapper->bp_comparative[s]       = (sliding_window) ? NULL : scs[s]->energy_bp;
            sc_wrapper->bp_local_comparative[s] =
              (sliding_window) ? scs[s]->energy_bp_local : NULL;
            sc_wrapper->user_cb_comparative[s]    = scs[s]->f;
            sc_wrapper->user_data_comparative[s]  = scs[s]->data;

            if (scs[s]->energy_up)
              provides_sc_up = 1;

            if (sliding_window) {
              if (scs[s]->energy_bp_local)
                provides_sc_bp = 1;
            } else if (scs[s]->energy_bp) {
              provides_sc_bp = 1;
            }

            if (scs[s]->f)
              provides_sc_user = 1;
          }
        }

        if (provides_sc_user) {
          sc_wrapper->pair_ext = &sc_hp_cb_ext_user_comparative;
          if (provides_sc_up) {
            sc_wrapper->pair_ext = &sc_hp_cb_ext_up_user_comparative;

            if (provides_sc_bp)
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_hp_cb_up_bp_local_user_comparative :
                                 &sc_hp_cb_up_bp_user_comparative;
            else
              sc_wrapper->pair = &sc_hp_cb_up_user_comparative;
          } else if (provides_sc_bp) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_hp_cb_bp_local_user_comparative :
                               &sc_hp_cb_bp_user_comparative;
          } else {
            sc_wrapper->pair = &sc_hp_cb_user_comparative;
          }
        } else {
          if (provides_sc_up) {
            sc_wrapper->pair_ext = &sc_hp_cb_ext_up_comparative;
            if (provides_sc_bp) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_hp_cb_up_bp_local_comparative :
                                 &sc_hp_cb_up_bp_comparative;
            } else {
              sc_wrapper->pair = &sc_hp_cb_up_comparative;
            }
          } else if (provides_sc_bp) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_hp_cb_bp_local_comparative :
                               &sc_hp_cb_bp_comparative;
          }
        }
      }

      break;
  }
}


PRIVATE INLINE void
free_sc_hp(struct sc_hp_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->bp_comparative);
  free(sc_wrapper->bp_local_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}
