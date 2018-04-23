struct sc_wrapper_ml;

typedef int (sc_pair)(int                   i,
                      int                   j,
                      struct sc_wrapper_ml  *data);

struct sc_wrapper_ml {
  int                     *idx;
  int                     **up;
  int                     *bp;

  sc_pair                 *pair;
  sc_pair                 *pair5;
  sc_pair                 *pair3;
  sc_pair                 *pair53;

  vrna_callback_sc_energy *user_cb;
  void                    *user_data;
};


PRIVATE INLINE int
sc_ml_pair_bp(int                   i,
              int                   j,
              struct sc_wrapper_ml  *data)
{
  return data->bp[data->idx[j] + i];
}


PRIVATE INLINE int
sc_ml_pair5_up(int                  i,
               int                  j,
               struct sc_wrapper_ml *data)
{
  return data->up[i + 1][1];
}


PRIVATE INLINE int
sc_ml_pair3_up(int                  i,
               int                  j,
               struct sc_wrapper_ml *data)
{
  return data->up[j - 1][1];
}


PRIVATE INLINE int
sc_ml_pair5_bp_up(int                   i,
                  int                   j,
                  struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair5_up(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_up(int                   i,
                  int                   j,
                  struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair3_up(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_up(int                  i,
                   int                  j,
                   struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair5_up(i, j, data) +
         sc_ml_pair3_up(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_up(int                   i,
                int                   j,
                struct sc_wrapper_ml  *data)
{
  return sc_ml_pair5_up(i, j, data) +
         sc_ml_pair3_up(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair_user(int                   i,
                int                   j,
                struct sc_wrapper_ml  *data)
{
  return data->user_cb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair_bp_user(int                  i,
                   int                  j,
                   struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         data->user_cb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair5_user(int                  i,
                 int                  j,
                 struct sc_wrapper_ml *data)
{
  return data->user_cb(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_user(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp(i, j, data) +
         data->user_cb(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair5_up_user(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair5_up(i, j, data) +
         data->user_cb(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_up_user(int                  i,
                       int                  j,
                       struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair5_up(i, j, data) +
         data->user_cb(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair3_user(int                  i,
                 int                  j,
                 struct sc_wrapper_ml *data)
{
  return data->user_cb(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_user(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp(i, j, data) +
         data->user_cb(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair3_up_user(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair3_up(i, j, data) +
         data->user_cb(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_up_user(int                  i,
                       int                  j,
                       struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair3_up(i, j, data) +
         data->user_cb(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair53_user(int                   i,
                  int                   j,
                  struct sc_wrapper_ml  *data)
{
  return data->user_cb(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_user(int                  i,
                     int                  j,
                     struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         data->user_cb(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair53_up_user(int                  i,
                     int                  j,
                     struct sc_wrapper_ml *data)
{
  return sc_ml_pair53_up(i, j, data) +
         data->user_cb(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_up_user(int                   i,
                        int                   j,
                        struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair53_up(i, j, data) +
         data->user_cb(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE void
init_sc_wrapper(vrna_fold_compound_t  *fc,
                struct sc_wrapper_ml  *sc_wrapper)
{
  vrna_sc_t *sc;

  sc = fc->sc;

  sc_wrapper->idx = fc->jindx;
  sc_wrapper->up  = NULL;
  sc_wrapper->bp  = NULL;

  sc_wrapper->pair    = NULL;
  sc_wrapper->pair5   = NULL;
  sc_wrapper->pair3   = NULL;
  sc_wrapper->pair53  = NULL;

  sc_wrapper->user_cb   = NULL;
  sc_wrapper->user_data = NULL;

  if (sc) {
    sc_wrapper->up        = sc->energy_up;
    sc_wrapper->bp        = sc->energy_bp;
    sc_wrapper->user_cb   = sc->f;
    sc_wrapper->user_data = sc->data;

    if (sc->f) {
      if (sc->energy_bp) {
        sc_wrapper->pair = &sc_ml_pair_bp_user;
        if (sc->energy_up) {
          sc_wrapper->pair5   = &sc_ml_pair5_bp_up_user;
          sc_wrapper->pair3   = &sc_ml_pair3_bp_up_user;
          sc_wrapper->pair53  = &sc_ml_pair53_bp_up_user;
        } else {
          sc_wrapper->pair5   = &sc_ml_pair5_bp_user;
          sc_wrapper->pair3   = &sc_ml_pair3_bp_user;
          sc_wrapper->pair53  = &sc_ml_pair53_bp_user;
        }
      } else if (sc->energy_up) {
        sc_wrapper->pair    = &sc_ml_pair_user;
        sc_wrapper->pair5   = &sc_ml_pair5_up_user;
        sc_wrapper->pair3   = &sc_ml_pair3_up_user;
        sc_wrapper->pair53  = &sc_ml_pair53_up_user;
      } else {
        sc_wrapper->pair    = &sc_ml_pair_user;
        sc_wrapper->pair5   = &sc_ml_pair5_user;
        sc_wrapper->pair3   = &sc_ml_pair3_user;
        sc_wrapper->pair53  = &sc_ml_pair53_user;
      }
    } else if (sc->energy_bp) {
      sc_wrapper->pair = &sc_ml_pair_bp;
      if (sc->energy_up) {
        sc_wrapper->pair5   = &sc_ml_pair5_bp_up;
        sc_wrapper->pair3   = &sc_ml_pair3_bp_up;
        sc_wrapper->pair53  = &sc_ml_pair53_bp_up;
      } else {
        sc_wrapper->pair5   = &sc_ml_pair_bp;
        sc_wrapper->pair3   = &sc_ml_pair_bp;
        sc_wrapper->pair53  = &sc_ml_pair_bp;
      }
    } else if (sc->energy_up) {
      sc_wrapper->pair5   = &sc_ml_pair5_up;
      sc_wrapper->pair3   = &sc_ml_pair3_up;
      sc_wrapper->pair53  = &sc_ml_pair53_up;
    }
  }
}


