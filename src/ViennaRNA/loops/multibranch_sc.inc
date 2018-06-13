struct sc_wrapper_ml;

typedef int (sc_pair)(int                   i,
                      int                   j,
                      struct sc_wrapper_ml  *data);

typedef int (sc_red)(int                  i,
                     int                  j,
                     int                  k,
                     int                  l,
                     struct sc_wrapper_ml *data);

struct sc_wrapper_ml {
  int                     n_seq;
  unsigned int            **a2s;

  int                     *idx;
  int                     **up;   /* unpaired constraints */
  int                     ***up_comparative;
  int                     *bp;    /* base pair constraints */
  int                     **bp_comparative;
  int                     **bp_local;
  int                     ***bp_local_comparative;
  int                     *stack; /* stack constraints (per nucleotide) */
  int                     **stack_comparative;

  sc_pair                 *pair;
  sc_pair                 *pair5;
  sc_pair                 *pair3;
  sc_pair                 *pair53;

  sc_red                  *red_stem;
  sc_red                  *red_ml;

  sc_red                  *decomp_ml;
  sc_red                  *coaxial_cls;
  sc_red                  *coaxial_enc;

  vrna_callback_sc_energy *user_cb;
  void                    *user_data;

  vrna_callback_sc_energy **user_cb_comparative;
  void                    **user_data_comparative;
};


PRIVATE INLINE int
sc_ml_pair_bp(int                   i,
              int                   j,
              struct sc_wrapper_ml  *data)
{
  return data->bp[data->idx[j] + i];
}


PRIVATE INLINE int
sc_ml_pair_bp_local(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return data->bp_local[i][j - i];
}


PRIVATE INLINE int
sc_ml_pair_bp_comparative(int                   i,
                          int                   j,
                          struct sc_wrapper_ml  *data)
{
  unsigned int  s;
  int           e, **sc_bp;

  e     = 0;
  sc_bp = data->bp_comparative;

  for (s = 0; s < data->n_seq; s++)
    if (sc_bp[s])
      e += sc_bp[s][data->idx[j] + i];

  return e;
}


PRIVATE INLINE int
sc_ml_pair_bp_local_comparative(int                   i,
                                int                   j,
                                struct sc_wrapper_ml  *data)
{
  unsigned int  s;
  int           e, ***sc_bp;

  e     = 0;
  sc_bp = data->bp_local_comparative;

  for (s = 0; s < data->n_seq; s++)
    if (sc_bp[s])
      e += sc_bp[s][i][j - i];

  return e;
}


PRIVATE INLINE int
sc_ml_pair5_up(int                  i,
               int                  j,
               struct sc_wrapper_ml *data)
{
  return data->up[i + 1][1];
}


PRIVATE INLINE int
sc_ml_pair5_up_comparative(int                  i,
                           int                  j,
                           struct sc_wrapper_ml *data)
{
  unsigned int  s, length, **a2s;
  int           e, ***sc_up;

  e     = 0;
  a2s   = data->a2s;
  sc_up = data->up_comparative;

  for (s = 0; s < data->n_seq; s++)
    if (sc_up[s]) {
      length  = a2s[s][i + 1] - a2s[s][i];
      e       += sc_up[s][a2s[s][i + 1]][length];
    }

  return e;
}


PRIVATE INLINE int
sc_ml_pair3_up(int                  i,
               int                  j,
               struct sc_wrapper_ml *data)
{
  return data->up[j - 1][1];
}


PRIVATE INLINE int
sc_ml_pair3_up_comparative(int                  i,
                           int                  j,
                           struct sc_wrapper_ml *data)
{
  unsigned int  s, length, **a2s;
  int           e, ***sc_up;

  e     = 0;
  a2s   = data->a2s;
  sc_up = data->up_comparative;

  for (s = 0; s < data->n_seq; s++)
    if (sc_up[s]) {
      length  = a2s[s][j] - a2s[s][j - 1];
      e       += sc_up[s][a2s[s][j - 1]][length];
    }

  return e;
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
sc_ml_pair5_bp_local_up(int                   i,
                        int                   j,
                        struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair5_up(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_up_comparative(int                   i,
                              int                   j,
                              struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair5_up_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_local_up_comparative(int                   i,
                                    int                   j,
                                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair5_up_comparative(i, j, data);
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
sc_ml_pair3_bp_local_up(int                   i,
                        int                   j,
                        struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair3_up(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_up_comparative(int                   i,
                              int                   j,
                              struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair3_up_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_local_up_comparative(int                   i,
                                    int                   j,
                                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair3_up_comparative(i, j, data);
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
sc_ml_pair53_bp_local_up(int                  i,
                         int                  j,
                         struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair5_up(i, j, data) +
         sc_ml_pair3_up(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_up_comparative(int                  i,
                               int                  j,
                               struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair5_up_comparative(i, j, data) +
         sc_ml_pair3_up_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_local_up_comparative(int                  i,
                                     int                  j,
                                     struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair5_up_comparative(i, j, data) +
         sc_ml_pair3_up_comparative(i, j, data);
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
sc_ml_pair53_up_comparative(int                   i,
                            int                   j,
                            struct sc_wrapper_ml  *data)
{
  return sc_ml_pair5_up_comparative(i, j, data) +
         sc_ml_pair3_up_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair_user(int                   i,
                int                   j,
                struct sc_wrapper_ml  *data)
{
  return data->user_cb(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair_user_comparative(int                   i,
                            int                   j,
                            struct sc_wrapper_ml  *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, i + 1, j - 1,
                                        VRNA_DECOMP_PAIR_ML,
                                        data->user_data_comparative[s]);

  return e;
}


PRIVATE INLINE int
sc_ml_pair_bp_user(int                  i,
                   int                  j,
                   struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair_bp_local_user(int                  i,
                         int                  j,
                         struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair_bp_user_comparative(int                  i,
                               int                  j,
                               struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair_bp_local_user_comparative(int                  i,
                                     int                  j,
                                     struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_user(int                  i,
                 int                  j,
                 struct sc_wrapper_ml *data)
{
  return data->user_cb(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair5_user_comparative(int                  i,
                             int                  j,
                             struct sc_wrapper_ml *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, i + 2, j - 1,
                                        VRNA_DECOMP_PAIR_ML,
                                        data->user_data);

  return e;
}


PRIVATE INLINE int
sc_ml_pair5_bp_user(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair5_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_local_user(int                   i,
                          int                   j,
                          struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair5_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_user_comparative(int                   i,
                                int                   j,
                                struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair5_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_local_user_comparative(int                   i,
                                      int                   j,
                                      struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair5_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_up_user(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair5_up(i, j, data) +
         sc_ml_pair5_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_up_user_comparative(int                   i,
                                int                   j,
                                struct sc_wrapper_ml  *data)
{
  return sc_ml_pair5_up_comparative(i, j, data) +
         sc_ml_pair5_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_up_user(int                  i,
                       int                  j,
                       struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair5_up(i, j, data) +
         sc_ml_pair5_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_local_up_user(int                  i,
                             int                  j,
                             struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair5_up(i, j, data) +
         sc_ml_pair5_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_up_user_comparative(int                  i,
                                   int                  j,
                                   struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair5_up_comparative(i, j, data) +
         sc_ml_pair5_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair5_bp_local_up_user_comparative(int                  i,
                                         int                  j,
                                         struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair5_up_comparative(i, j, data) +
         sc_ml_pair5_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_user(int                  i,
                 int                  j,
                 struct sc_wrapper_ml *data)
{
  return data->user_cb(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair3_user_comparative(int                  i,
                             int                  j,
                             struct sc_wrapper_ml *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, i + 1, j - 2,
                                        VRNA_DECOMP_PAIR_ML,
                                        data->user_data);

  return e;
}


PRIVATE INLINE int
sc_ml_pair3_bp_user(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair3_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_local_user(int                   i,
                          int                   j,
                          struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair3_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_user_comparative(int                   i,
                                int                   j,
                                struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair3_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_local_user_comparative(int                   i,
                                      int                   j,
                                      struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair3_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_up_user(int                   i,
                    int                   j,
                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair3_up(i, j, data) +
         sc_ml_pair3_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_up_user_comparative(int                   i,
                                int                   j,
                                struct sc_wrapper_ml  *data)
{
  return sc_ml_pair3_up_comparative(i, j, data) +
         sc_ml_pair3_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_up_user(int                  i,
                       int                  j,
                       struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair3_up(i, j, data) +
         sc_ml_pair3_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_local_up_user(int                  i,
                             int                  j,
                             struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair3_up(i, j, data) +
         sc_ml_pair3_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_up_user_comparative(int                  i,
                                   int                  j,
                                   struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair3_up_comparative(i, j, data) +
         sc_ml_pair3_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair3_bp_local_up_user_comparative(int                  i,
                                         int                  j,
                                         struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair3_up_comparative(i, j, data) +
         sc_ml_pair3_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_user(int                   i,
                  int                   j,
                  struct sc_wrapper_ml  *data)
{
  return data->user_cb(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, data->user_data);
}


PRIVATE INLINE int
sc_ml_pair53_user_comparative(int                   i,
                              int                   j,
                              struct sc_wrapper_ml  *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, i + 2, j - 2,
                                        VRNA_DECOMP_PAIR_ML,
                                        data->user_data);

  return e;
}


PRIVATE INLINE int
sc_ml_pair53_bp_user(int                  i,
                     int                  j,
                     struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair53_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_local_user(int                  i,
                           int                  j,
                           struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair53_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_user_comparative(int                  i,
                                 int                  j,
                                 struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair53_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_local_user_comparative(int                  i,
                                       int                  j,
                                       struct sc_wrapper_ml *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair53_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_up_user(int                  i,
                     int                  j,
                     struct sc_wrapper_ml *data)
{
  return sc_ml_pair53_up(i, j, data) +
         sc_ml_pair53_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_up_user_comparative(int                  i,
                                 int                  j,
                                 struct sc_wrapper_ml *data)
{
  return sc_ml_pair53_up_comparative(i, j, data) +
         sc_ml_pair53_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_up_user(int                   i,
                        int                   j,
                        struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp(i, j, data) +
         sc_ml_pair53_up(i, j, data) +
         sc_ml_pair53_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_local_up_user(int                   i,
                              int                   j,
                              struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local(i, j, data) +
         sc_ml_pair53_up(i, j, data) +
         sc_ml_pair53_user(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_up_user_comparative(int                   i,
                                    int                   j,
                                    struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_comparative(i, j, data) +
         sc_ml_pair53_up_comparative(i, j, data) +
         sc_ml_pair53_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_pair53_bp_local_up_user_comparative(int                   i,
                                          int                   j,
                                          struct sc_wrapper_ml  *data)
{
  return sc_ml_pair_bp_local_comparative(i, j, data) +
         sc_ml_pair53_up_comparative(i, j, data) +
         sc_ml_pair53_user_comparative(i, j, data);
}


PRIVATE INLINE int
sc_ml_reduce_stem_user(int                  i,
                       int                  j,
                       int                  k,
                       int                  l,
                       struct sc_wrapper_ml *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_STEM,
                       data->user_data);
}


PRIVATE INLINE int
sc_ml_reduce_stem_user_comparative(int                  i,
                                   int                  j,
                                   int                  k,
                                   int                  l,
                                   struct sc_wrapper_ml *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, k, l,
                                        VRNA_DECOMP_ML_STEM,
                                        data->user_data_comparative[s]);

  return e;
}


PRIVATE INLINE int
sc_ml_reduce_ml(int                   i,
                int                   j,
                int                   k,
                int                   l,
                struct sc_wrapper_ml  *data)
{
  unsigned int  u1, u2;
  int           e;

  e   = 0;
  u1  = k - i;
  u2  = j - l;

  if (u1)
    e += data->up[i][u1];

  if (u2)
    e += data->up[l + 1][u2];

  return e;
}


PRIVATE INLINE int
sc_ml_reduce_ml_comparative(int                   i,
                            int                   j,
                            int                   k,
                            int                   l,
                            struct sc_wrapper_ml  *data)
{
  unsigned int  s, u1, u2, **a2s;
  int           e;

  e   = 0;
  a2s = data->a2s;

  for (s = 0; s < data->n_seq; s++) {
    u1  = a2s[s][k] - a2s[s][i];
    u2  = a2s[s][j] - a2s[s][l];

    if (u1)
      e += data->up[a2s[s][i]][u1];

    if (u2)
      e += data->up[a2s[s][l] + 1][u2];
  }
  return e;
}


PRIVATE INLINE int
sc_ml_reduce_ml_user(int                  i,
                     int                  j,
                     int                  k,
                     int                  l,
                     struct sc_wrapper_ml *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_ML,
                       data->user_data);
}


PRIVATE INLINE int
sc_ml_reduce_ml_user_comparative(int                  i,
                                 int                  j,
                                 int                  k,
                                 int                  l,
                                 struct sc_wrapper_ml *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, k, l,
                                        VRNA_DECOMP_ML_ML,
                                        data->user_data_comparative[s]);

  return e;
}


PRIVATE INLINE int
sc_ml_reduce_ml_up_user(int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        struct sc_wrapper_ml  *data)
{
  return sc_ml_reduce_ml(i, j, k, l, data) +
         sc_ml_reduce_ml_user(i, j, k, l, data);
}


PRIVATE INLINE int
sc_ml_reduce_ml_up_user_comparative(int                   i,
                                    int                   j,
                                    int                   k,
                                    int                   l,
                                    struct sc_wrapper_ml  *data)
{
  return sc_ml_reduce_ml_comparative(i, j, k, l, data) +
         sc_ml_reduce_ml_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE int
sc_ml_decomp_ml_user(int                  i,
                     int                  j,
                     int                  k,
                     int                  l,
                     struct sc_wrapper_ml *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_ML_ML,
                       data->user_data);
}


PRIVATE INLINE int
sc_ml_decomp_ml_user_comparative(int                  i,
                                 int                  j,
                                 int                  k,
                                 int                  l,
                                 struct sc_wrapper_ml *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, k, l,
                                        VRNA_DECOMP_ML_ML_ML,
                                        data->user_data_comparative[s]);

  return e;
}


PRIVATE INLINE int
sc_ml_coax_stack(int                  i,
                 int                  j,
                 int                  k,
                 int                  l,
                 struct sc_wrapper_ml *data)
{
  return data->stack[i] +
         data->stack[j] +
         data->stack[k] +
         data->stack[l];
}


PRIVATE INLINE int
sc_ml_coax_stack_comparative(int                  i,
                             int                  j,
                             int                  k,
                             int                  l,
                             struct sc_wrapper_ml *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->stack_comparative[s]) {
      e += data->stack_comparative[s][data->a2s[s][i]] +
           data->stack_comparative[s][data->a2s[s][j]] +
           data->stack_comparative[s][data->a2s[s][k]] +
           data->stack_comparative[s][data->a2s[s][l]];
    }

  return e;
}


PRIVATE INLINE int
sc_ml_coax_closing_user(int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        struct sc_wrapper_ml  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_COAXIAL,
                       data->user_data);
}


PRIVATE INLINE int
sc_ml_coax_closing_user_comparative(int                   i,
                                    int                   j,
                                    int                   k,
                                    int                   l,
                                    struct sc_wrapper_ml  *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, k, l,
                                        VRNA_DECOMP_ML_COAXIAL,
                                        data->user_data_comparative[s]);

  return e;
}


PRIVATE INLINE int
sc_ml_coax_closing_stack_user(int                   i,
                              int                   j,
                              int                   k,
                              int                   l,
                              struct sc_wrapper_ml  *data)
{
  return sc_ml_coax_stack(i, j, k, l, data) +
         sc_ml_coax_closing_user(i, j, k, l, data);
}


PRIVATE INLINE int
sc_ml_coax_closing_stack_user_comparative(int                   i,
                                          int                   j,
                                          int                   k,
                                          int                   l,
                                          struct sc_wrapper_ml  *data)
{
  return sc_ml_coax_stack_comparative(i, j, k, l, data) +
         sc_ml_coax_closing_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE int
sc_ml_coax_enclosed_user(int                  i,
                         int                  j,
                         int                  k,
                         int                  l,
                         struct sc_wrapper_ml *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_ML_COAXIAL_ENC,
                       data->user_data);
}


PRIVATE INLINE int
sc_ml_coax_enclosed_user_comparative(int                  i,
                                     int                  j,
                                     int                  k,
                                     int                  l,
                                     struct sc_wrapper_ml *data)
{
  unsigned int  s;
  int           e;

  e = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e += data->user_cb_comparative[s](i, j, k, l,
                                        VRNA_DECOMP_ML_COAXIAL_ENC,
                                        data->user_data_comparative[s]);

  return e;
}


PRIVATE INLINE int
sc_ml_coax_enclosed_stack_user(int                  i,
                               int                  j,
                               int                  k,
                               int                  l,
                               struct sc_wrapper_ml *data)
{
  return sc_ml_coax_stack(i, j, k, l, data) +
         sc_ml_coax_enclosed_user(i, j, k, l, data);
}


PRIVATE INLINE int
sc_ml_coax_enclosed_stack_user_comparative(int                  i,
                                           int                  j,
                                           int                  k,
                                           int                  l,
                                           struct sc_wrapper_ml *data)
{
  return sc_ml_coax_stack_comparative(i, j, k, l, data) +
         sc_ml_coax_enclosed_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE void
init_sc_wrapper(vrna_fold_compound_t  *fc,
                struct sc_wrapper_ml  *sc_wrapper)
{
  unsigned char sliding_window;
  vrna_sc_t     *sc, **scs;

  sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;

  sc_wrapper->idx                   = fc->jindx;
  sc_wrapper->n_seq                 = 1;
  sc_wrapper->a2s                   = NULL;
  sc_wrapper->up                    = NULL;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->bp                    = NULL;
  sc_wrapper->bp_comparative        = NULL;
  sc_wrapper->bp_local              = NULL;
  sc_wrapper->bp_local_comparative  = NULL;
  sc_wrapper->stack                 = NULL;
  sc_wrapper->stack_comparative     = NULL;

  sc_wrapper->pair    = NULL;
  sc_wrapper->pair5   = NULL;
  sc_wrapper->pair3   = NULL;
  sc_wrapper->pair53  = NULL;

  sc_wrapper->red_stem  = NULL;
  sc_wrapper->red_ml    = NULL;

  sc_wrapper->decomp_ml   = NULL;
  sc_wrapper->coaxial_cls = NULL;
  sc_wrapper->coaxial_enc = NULL;

  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data_comparative = NULL;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        unsigned int provides_sc_up, provides_sc_bp, provides_sc_stack, provides_sc_user;

        provides_sc_up    = 0;
        provides_sc_bp    = 0;
        provides_sc_stack = 0;
        provides_sc_user  = 0;

        sc_wrapper->up        = sc->energy_up;
        sc_wrapper->stack     = sc->energy_stack;
        sc_wrapper->user_cb   = sc->f;
        sc_wrapper->user_data = sc->data;

        if (fc->hc->type == VRNA_HC_WINDOW) {
          sc_wrapper->bp_local = sc->energy_bp_local;

          if (sc->energy_bp_local)
            provides_sc_bp = 1;
        } else {
          sc_wrapper->bp = sc->energy_bp;

          if (sc->energy_bp)
            provides_sc_bp = 1;
        }

        if (sc->energy_up)
          provides_sc_up = 1;

        if (sc->energy_stack)
          provides_sc_stack = 1;

        if (sc->f)
          provides_sc_user = 1;

        if (provides_sc_user) {
          sc_wrapper->red_stem  = &sc_ml_reduce_stem_user;
          sc_wrapper->red_ml    = &sc_ml_reduce_ml_user;
          sc_wrapper->decomp_ml = &sc_ml_decomp_ml_user;

          if (provides_sc_stack) {
            sc_wrapper->coaxial_cls = &sc_ml_coax_closing_stack_user;
            sc_wrapper->coaxial_enc = &sc_ml_coax_enclosed_stack_user;
          } else {
            sc_wrapper->coaxial_cls = &sc_ml_coax_closing_user;
            sc_wrapper->coaxial_enc = &sc_ml_coax_enclosed_user;
          }

          if (provides_sc_bp) {
            if (sliding_window) {
              sc_wrapper->pair = &sc_ml_pair_bp_local_user;
              if (provides_sc_up) {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_local_up_user;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_local_up_user;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_local_up_user;

                sc_wrapper->red_ml = &sc_ml_reduce_ml_up_user;
              } else {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_local_user;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_local_user;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_local_user;
              }
            } else {
              sc_wrapper->pair = &sc_ml_pair_bp_user;
              if (provides_sc_up) {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_up_user;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_up_user;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_up_user;

                sc_wrapper->red_ml = &sc_ml_reduce_ml_up_user;
              } else {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_user;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_user;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_user;
              }
            }
          } else if (provides_sc_up) {
            sc_wrapper->pair    = &sc_ml_pair_user;
            sc_wrapper->pair5   = &sc_ml_pair5_up_user;
            sc_wrapper->pair3   = &sc_ml_pair3_up_user;
            sc_wrapper->pair53  = &sc_ml_pair53_up_user;

            sc_wrapper->red_ml = &sc_ml_reduce_ml_up_user;
          } else {
            sc_wrapper->pair    = &sc_ml_pair_user;
            sc_wrapper->pair5   = &sc_ml_pair5_user;
            sc_wrapper->pair3   = &sc_ml_pair3_user;
            sc_wrapper->pair53  = &sc_ml_pair53_user;
          }
        } else {
          if (provides_sc_stack) {
            sc_wrapper->coaxial_cls = &sc_ml_coax_stack;
            sc_wrapper->coaxial_enc = &sc_ml_coax_stack;
          }

          if (provides_sc_bp) {
            if (sliding_window) {
              sc_wrapper->pair = &sc_ml_pair_bp_local;
              if (provides_sc_up) {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_local_up;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_local_up;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_local_up;

                sc_wrapper->red_ml = &sc_ml_reduce_ml;
              } else {
                sc_wrapper->pair5   = &sc_ml_pair_bp_local;
                sc_wrapper->pair3   = &sc_ml_pair_bp_local;
                sc_wrapper->pair53  = &sc_ml_pair_bp_local;
              }
            } else {
              sc_wrapper->pair = &sc_ml_pair_bp;
              if (provides_sc_up) {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_up;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_up;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_up;

                sc_wrapper->red_ml = &sc_ml_reduce_ml;
              } else {
                sc_wrapper->pair5   = &sc_ml_pair_bp;
                sc_wrapper->pair3   = &sc_ml_pair_bp;
                sc_wrapper->pair53  = &sc_ml_pair_bp;
              }
            }
          } else if (provides_sc_up) {
            sc_wrapper->pair5   = &sc_ml_pair5_up;
            sc_wrapper->pair3   = &sc_ml_pair3_up;
            sc_wrapper->pair53  = &sc_ml_pair53_up;

            sc_wrapper->red_ml = &sc_ml_reduce_ml;
          }
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      scs               = fc->scs;
      sc_wrapper->n_seq = fc->n_seq;
      sc_wrapper->a2s   = fc->a2s;

      if (scs) {
        unsigned int s, provides_sc_up, provides_sc_bp, provides_sc_user_cb, provides_sc_stack;

        provides_sc_up      = 0;
        provides_sc_bp      = 0;
        provides_sc_user_cb = 0;
        provides_sc_stack   = 0;

        sc_wrapper->up_comparative = (int ***)vrna_alloc(
          sizeof(int **) *
          fc->n_seq);
        sc_wrapper->bp_comparative = (int **)vrna_alloc(
          sizeof(int *) *
          fc->n_seq);
        sc_wrapper->bp_local_comparative = (int ***)vrna_alloc(
          sizeof(int **) *
          fc->n_seq);
        sc_wrapper->stack_comparative = (int **)vrna_alloc(
          sizeof(int *) *
          fc->n_seq);
        sc_wrapper->user_cb_comparative = (vrna_callback_sc_energy **)vrna_alloc(
          sizeof(vrna_callback_sc_energy *) *
          fc->n_seq);
        sc_wrapper->user_data_comparative = (void **)vrna_alloc(
          sizeof(void *) *
          fc->n_seq);

        for (s = 0; s < fc->n_seq; s++) {
          if (scs[s]) {
            sc_wrapper->up_comparative[s]       = scs[s]->energy_up;
            sc_wrapper->bp_comparative[s]       = (sliding_window) ? NULL : scs[s]->energy_bp;
            sc_wrapper->bp_local_comparative[s] =
              (sliding_window) ? scs[s]->energy_bp_local : NULL;
            sc_wrapper->stack_comparative[s]      = scs[s]->energy_stack;
            sc_wrapper->user_cb_comparative[s]    = scs[s]->f;
            sc_wrapper->user_data_comparative[s]  = scs[s]->data;

            if (scs[s]->energy_up)
              provides_sc_up = 1;

            if (((sliding_window) && (scs[s]->energy_bp_local)) || (scs[s]->energy_bp))
              provides_sc_bp = 1;

            if (scs[s]->f)
              provides_sc_user_cb = 1;

            if (scs[s]->energy_stack)
              provides_sc_stack = 1;
          }
        }

        /* bind the correct wrapper functions */
        if (provides_sc_user_cb) {
          sc_wrapper->red_stem  = &sc_ml_reduce_stem_user_comparative;
          sc_wrapper->red_ml    = &sc_ml_reduce_ml_user_comparative;
          sc_wrapper->decomp_ml = &sc_ml_decomp_ml_user_comparative;

          if (provides_sc_stack) {
            sc_wrapper->coaxial_cls = &sc_ml_coax_closing_stack_user_comparative;
            sc_wrapper->coaxial_enc = &sc_ml_coax_enclosed_stack_user_comparative;
          } else {
            sc_wrapper->coaxial_cls = &sc_ml_coax_closing_user_comparative;
            sc_wrapper->coaxial_enc = &sc_ml_coax_enclosed_user_comparative;
          }

          if (provides_sc_bp) {
            if (sliding_window) {
              sc_wrapper->pair = &sc_ml_pair_bp_local_user_comparative;
              if (provides_sc_up) {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_local_up_user_comparative;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_local_up_user_comparative;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_local_up_user_comparative;

                sc_wrapper->red_ml = &sc_ml_reduce_ml_up_user_comparative;
              } else {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_local_user_comparative;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_local_user_comparative;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_local_user_comparative;
              }
            } else {
              sc_wrapper->pair = &sc_ml_pair_bp_user_comparative;
              if (provides_sc_up) {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_up_user_comparative;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_up_user_comparative;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_up_user_comparative;

                sc_wrapper->red_ml = &sc_ml_reduce_ml_up_user_comparative;
              } else {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_user_comparative;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_user_comparative;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_user_comparative;
              }
            }
          } else if (provides_sc_up) {
            sc_wrapper->pair    = &sc_ml_pair_user_comparative;
            sc_wrapper->pair5   = &sc_ml_pair5_up_user_comparative;
            sc_wrapper->pair3   = &sc_ml_pair3_up_user_comparative;
            sc_wrapper->pair53  = &sc_ml_pair53_up_user_comparative;

            sc_wrapper->red_ml = &sc_ml_reduce_ml_up_user_comparative;
          } else {
            sc_wrapper->pair    = &sc_ml_pair_user_comparative;
            sc_wrapper->pair5   = &sc_ml_pair5_user_comparative;
            sc_wrapper->pair3   = &sc_ml_pair3_user_comparative;
            sc_wrapper->pair53  = &sc_ml_pair53_user_comparative;
          }
        } else {
          if (provides_sc_stack) {
            sc_wrapper->coaxial_cls = &sc_ml_coax_stack_comparative;
            sc_wrapper->coaxial_enc = &sc_ml_coax_stack_comparative;
          }

          if (provides_sc_bp) {
            if (sliding_window) {
              sc_wrapper->pair = &sc_ml_pair_bp_local_comparative;

              if (provides_sc_up) {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_local_up_comparative;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_local_up_comparative;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_local_up_comparative;

                sc_wrapper->red_ml = &sc_ml_reduce_ml_comparative;
              } else {
                sc_wrapper->pair5   = &sc_ml_pair_bp_local_comparative;
                sc_wrapper->pair3   = &sc_ml_pair_bp_local_comparative;
                sc_wrapper->pair53  = &sc_ml_pair_bp_local_comparative;
              }
            } else {
              sc_wrapper->pair = &sc_ml_pair_bp_comparative;

              if (provides_sc_up) {
                sc_wrapper->pair5   = &sc_ml_pair5_bp_up_comparative;
                sc_wrapper->pair3   = &sc_ml_pair3_bp_up_comparative;
                sc_wrapper->pair53  = &sc_ml_pair53_bp_up_comparative;

                sc_wrapper->red_ml = &sc_ml_reduce_ml_comparative;
              } else {
                sc_wrapper->pair5   = &sc_ml_pair_bp_comparative;
                sc_wrapper->pair3   = &sc_ml_pair_bp_comparative;
                sc_wrapper->pair53  = &sc_ml_pair_bp_comparative;
              }
            }
          } else if (provides_sc_up) {
            sc_wrapper->pair5   = &sc_ml_pair5_up_comparative;
            sc_wrapper->pair3   = &sc_ml_pair3_up_comparative;
            sc_wrapper->pair53  = &sc_ml_pair53_up_comparative;
            sc_wrapper->red_ml  = &sc_ml_reduce_ml_comparative;
          }
        }
      }

      break;
  }
}


PRIVATE INLINE void
free_sc_wrapper(struct sc_wrapper_ml *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->bp_comparative);
  free(sc_wrapper->stack_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}
