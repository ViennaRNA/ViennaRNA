#include <ViennaRNA/utils/basic.h>


struct sc_int_exp_dat;

typedef FLT_OR_DBL (*sc_int_exp_cb)(int                    i,
                                   int                    j,
                                   int                    k,
                                   int                    l,
                                   struct sc_int_exp_dat  *data);

struct sc_int_exp_dat {
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
  FLT_OR_DBL                  *stack;
  FLT_OR_DBL                  **stack_comparative;

  vrna_sc_exp_f user_cb;
  void                        *user_data;

  vrna_sc_exp_f *user_cb_comparative;
  void                        **user_data_comparative;

  sc_int_exp_cb               pair;
  sc_int_exp_cb               pair_ext;
};


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up(int                    i,
                 int                    j,
                 int                    k,
                 int                    l,
                 struct sc_int_exp_dat  *data)
{
  int         u1, u2;
  FLT_OR_DBL  sc;

  u1  = k - i - 1;
  u2  = j - l - 1;

  sc = 1.;

  if (u1 > 0)
    sc *= data->up[i + 1][u1];

  if (u2 > 0)
    sc *= data->up[l + 1][u2];

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_comparative(int                    i,
                             int                    j,
                             int                    k,
                             int                    l,
                             struct sc_int_exp_dat  *data)
{
  unsigned int  s, u1, u2;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->up_comparative[s]) {
      u1  = data->a2s[s][k - 1] - data->a2s[s][i];
      u2  = data->a2s[s][j - 1] - data->a2s[s][l];

      if (u1 != 0)
        sc *= data->up_comparative[s][data->a2s[s][i] + 1][u1];

      if (u2 != 0)
        sc *= data->up_comparative[s][data->a2s[s][l] + 1][u2];
    }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp(int                    i,
                 int                    j,
                 int                    k VRNA_UNUSED,
                 int                    l VRNA_UNUSED,
                 struct sc_int_exp_dat  *data)
{
  return data->bp[data->idx[j] + i];
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_comparative(int                    i,
                             int                    j,
                             int                    k VRNA_UNUSED,
                             int                    l VRNA_UNUSED,
                             struct sc_int_exp_dat  *data)
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
sc_int_exp_cb_bp_local(int                    i,
                       int                    j,
                       int                    k VRNA_UNUSED,
                       int                    l VRNA_UNUSED,
                       struct sc_int_exp_dat  *data)
{
  return data->bp_local[i][j - i];
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_comparative(int                    i,
                                   int                    j,
                                   int                    k VRNA_UNUSED,
                                   int                    l VRNA_UNUSED,
                                   struct sc_int_exp_dat  *data)
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
sc_int_exp_cb_stack(int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    struct sc_int_exp_dat *data)
{
  FLT_OR_DBL sc;

  sc = 1.;

  if ((i + 1 == k) && (l + 1 == j)) {
    sc *= data->stack[i] *
          data->stack[k] *
          data->stack[l] *
          data->stack[j];
  }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_stack_comparative(int                   i,
                                int                   j,
                                int                   k,
                                int                   l,
                                struct sc_int_exp_dat *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->stack_comparative[s]) {
      if ((data->a2s[s][k] == data->a2s[s][i] + 1) &&
          (data->a2s[s][j] == data->a2s[s][l] + 1)) {
        sc *= data->stack_comparative[s][data->a2s[s][i]] *
              data->stack_comparative[s][data->a2s[s][k]] *
              data->stack_comparative[s][data->a2s[s][l]] *
              data->stack_comparative[s][data->a2s[s][j]];
      }
    }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_user(int                    i,
                   int                    j,
                   int                    k,
                   int                    l,
                   struct sc_int_exp_dat  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_PAIR_IL,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_user_comparative(int                    i,
                               int                    j,
                               int                    k,
                               int                    l,
                               struct sc_int_exp_dat  *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s]) {
      sc *= data->user_cb_comparative[s](i, j, k, l,
                                         VRNA_DECOMP_PAIR_IL,
                                         data->user_data_comparative[s]);
    }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp(int                   i,
                    int                   j,
                    int                   k,
                    int                   l,
                    struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_comparative(int                   i,
                                int                   j,
                                int                   k,
                                int                   l,
                                struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_bp_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local(int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp_local(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_comparative(int                   i,
                                      int                   j,
                                      int                   k,
                                      int                   l,
                                      struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_bp_local_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_stack(int                    i,
                       int                    j,
                       int                    k,
                       int                    l,
                       struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_stack_comparative(int                    i,
                                   int                    j,
                                   int                    k,
                                   int                    l,
                                   struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_user(int                   i,
                      int                   j,
                      int                   k,
                      int                   l,
                      struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_user_comparative(int                   i,
                                  int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_stack(int                    i,
                       int                    j,
                       int                    k,
                       int                    l,
                       struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_stack_comparative(int                    i,
                                   int                    j,
                                   int                    k,
                                   int                    l,
                                   struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_stack(int                    i,
                             int                    j,
                             int                    k,
                             int                    l,
                             struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_stack_comparative(int                    i,
                                         int                    j,
                                         int                    k,
                                         int                    l,
                                         struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_user(int                   i,
                      int                   j,
                      int                   k,
                      int                   l,
                      struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_user_comparative(int                   i,
                                  int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_user(int                   i,
                            int                   j,
                            int                   k,
                            int                   l,
                            struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_user_comparative(int                   i,
                                        int                   j,
                                        int                   k,
                                        int                   l,
                                        struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_stack_user(int                    i,
                         int                    j,
                         int                    k,
                         int                    l,
                         struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_stack_user_comparative(int                    i,
                                     int                    j,
                                     int                    k,
                                     int                    l,
                                     struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_stack(int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_stack_comparative(int                   i,
                                      int                   j,
                                      int                   k,
                                      int                   l,
                                      struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_stack(int                   i,
                                int                   j,
                                int                   k,
                                int                   l,
                                struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_stack_comparative(int                   i,
                                            int                   j,
                                            int                   k,
                                            int                   l,
                                            struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_user(int                    i,
                         int                    j,
                         int                    k,
                         int                    l,
                         struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_user_comparative(int                    i,
                                     int                    j,
                                     int                    k,
                                     int                    l,
                                     struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_user(int                    i,
                               int                    j,
                               int                    k,
                               int                    l,
                               struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_user_comparative(int                    i,
                                           int                    j,
                                           int                    k,
                                           int                    l,
                                           struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_stack_user(int                   i,
                            int                   j,
                            int                   k,
                            int                   l,
                            struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_stack_user_comparative(int                   i,
                                        int                   j,
                                        int                   k,
                                        int                   l,
                                        struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_stack_user(int                   i,
                            int                   j,
                            int                   k,
                            int                   l,
                            struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_stack_user_comparative(int                   i,
                                        int                   j,
                                        int                   k,
                                        int                   l,
                                        struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_stack_user(int                   i,
                                  int                   j,
                                  int                   k,
                                  int                   l,
                                  struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_bp_local_stack_user_comparative(int                   i,
                                              int                   j,
                                              int                   k,
                                              int                   l,
                                              struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_stack_user(int                    i,
                               int                    j,
                               int                    k,
                               int                    l,
                               struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_stack_user_comparative(int                    i,
                                           int                    j,
                                           int                    k,
                                           int                    l,
                                           struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_bp_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_stack_user(int                    i,
                                     int                    j,
                                     int                    k,
                                     int                    l,
                                     struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up(i, j, k, l, data) *
         sc_int_exp_cb_bp_local(i, j, k, l, data) *
         sc_int_exp_cb_stack(i, j, k, l, data) *
         sc_int_exp_cb_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_up_bp_local_stack_user_comparative(int                    i,
                                                 int                    j,
                                                 int                    k,
                                                 int                    l,
                                                 struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_bp_local_comparative(i, j, k, l, data) *
         sc_int_exp_cb_stack_comparative(i, j, k, l, data) *
         sc_int_exp_cb_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up(int                    i,
                     int                    j,
                     int                    k,
                     int                    l,
                     struct sc_int_exp_dat  *data)
{
  int         u1, u2, u3;
  FLT_OR_DBL  sc;

  sc  = 1.;
  u1  = i - 1;
  u2  = k - j - 1;
  u3  = data->n - l;

  if (u1 > 0)
    sc *= data->up[1][u1];

  if (u2 > 0)
    sc *= data->up[j + 1][u2];

  if (u3 > 0)
    sc *= data->up[l + 1][u3];

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_comparative(int                    i,
                                 int                    j,
                                 int                    k,
                                 int                    l,
                                 struct sc_int_exp_dat  *data)
{
  unsigned int  s;
  int           u1, u2, u3;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->up_comparative[s]) {
      u1  = (i > 0) ? data->a2s[s][i - 1] - data->a2s[s][1] : 0;
      u2  = data->a2s[s][k - 1] - data->a2s[s][j];
      u3  = data->a2s[s][data->n] - data->a2s[s][l];

      if (u1 != 0)
        sc *= data->up_comparative[s][data->a2s[s][1]][u1];

      if (u2 != 0)
        sc *= data->up_comparative[s][data->a2s[s][j] + 1][u2];

      if (u3 != 0)
        sc *= data->up_comparative[s][data->a2s[s][l] + 1][u3];
    }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_stack(int                   i,
                        int                   j,
                        int                   k,
                        int                   l,
                        struct sc_int_exp_dat *data)
{
  FLT_OR_DBL sc;

  sc = 1.;

  if ((i == 1) && (j + 1 == k) && (l == (int)data->n)) {
    sc *= data->stack[i] *
          data->stack[k] *
          data->stack[l] *
          data->stack[j];
  }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_stack_comparative(int                   i,
                                    int                   j,
                                    int                   k,
                                    int                   l,
                                    struct sc_int_exp_dat *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->stack_comparative[s]) {
      if ((data->a2s[s][i] == 1) &&
          (data->a2s[s][j] + 1 == data->a2s[s][k]) &&
          (data->a2s[s][l] == data->a2s[s][data->n])) {
        sc *= data->stack_comparative[s][data->a2s[s][i]] *
              data->stack_comparative[s][data->a2s[s][k]] *
              data->stack_comparative[s][data->a2s[s][l]] *
              data->stack_comparative[s][data->a2s[s][j]];
      }
    }

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_user(int                    i,
                       int                    j,
                       int                    k,
                       int                    l,
                       struct sc_int_exp_dat  *data)
{
  return data->user_cb(i, j, k, l,
                       VRNA_DECOMP_PAIR_IL,
                       data->user_data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_user_comparative(int                    i,
                                   int                    j,
                                   int                    k,
                                   int                    l,
                                   struct sc_int_exp_dat  *data)
{
  unsigned int  s;
  FLT_OR_DBL    sc;

  sc = 1.;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      sc *= data->user_cb_comparative[s](i, j, k, l,
                                         VRNA_DECOMP_PAIR_IL,
                                         data->user_data_comparative[s]);

  return sc;
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_stack(int                    i,
                           int                    j,
                           int                    k,
                           int                    l,
                           struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_ext_up(i, j, k, l, data) *
         sc_int_exp_cb_ext_stack(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_stack_comparative(int                    i,
                                       int                    j,
                                       int                    k,
                                       int                    l,
                                       struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_ext_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_ext_stack_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_user(int                   i,
                          int                   j,
                          int                   k,
                          int                   l,
                          struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_ext_up(i, j, k, l, data) *
         sc_int_exp_cb_ext_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_user_comparative(int                   i,
                                      int                   j,
                                      int                   k,
                                      int                   l,
                                      struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_ext_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_ext_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_stack_user(int                    i,
                             int                    j,
                             int                    k,
                             int                    l,
                             struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_ext_stack(i, j, k, l, data) *
         sc_int_exp_cb_ext_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_stack_user_comparative(int                    i,
                                         int                    j,
                                         int                    k,
                                         int                    l,
                                         struct sc_int_exp_dat  *data)
{
  return sc_int_exp_cb_ext_stack_comparative(i, j, k, l, data) *
         sc_int_exp_cb_ext_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_stack_user(int                   i,
                                int                   j,
                                int                   k,
                                int                   l,
                                struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_ext_up(i, j, k, l, data) *
         sc_int_exp_cb_ext_stack(i, j, k, l, data) *
         sc_int_exp_cb_ext_user(i, j, k, l, data);
}


PRIVATE INLINE FLT_OR_DBL
sc_int_exp_cb_ext_up_stack_user_comparative(int                   i,
                                            int                   j,
                                            int                   k,
                                            int                   l,
                                            struct sc_int_exp_dat *data)
{
  return sc_int_exp_cb_ext_up_comparative(i, j, k, l, data) *
         sc_int_exp_cb_ext_stack_comparative(i, j, k, l, data) *
         sc_int_exp_cb_ext_user_comparative(i, j, k, l, data);
}


PRIVATE INLINE void
init_sc_int_exp(vrna_fold_compound_t  *fc,
                struct sc_int_exp_dat *sc_wrapper)
{
  unsigned char sliding_window;
  unsigned int  s, provides_sc_up, provides_sc_bp, provides_sc_stack, provides_sc_user;
  vrna_sc_t     *sc, **scs;

  if (fc->exp_matrices)
    sliding_window = (fc->exp_matrices->type == VRNA_MX_WINDOW) ? 1 : 0;
  else if ((fc->type == VRNA_FC_TYPE_SINGLE) && (fc->sc))
    sliding_window = (fc->sc->type == VRNA_SC_WINDOW) ? 1 : 0;
  else if (fc->hc)
    sliding_window = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  else
    sliding_window = 0;

  provides_sc_up    = 0;
  provides_sc_bp    = 0;
  provides_sc_stack = 0;
  provides_sc_user  = 0;

  sc_wrapper->n                     = fc->length;
  sc_wrapper->n_seq                 = 1;
  sc_wrapper->a2s                   = NULL;
  sc_wrapper->idx                   = fc->jindx;
  sc_wrapper->up                    = NULL;
  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->bp                    = NULL;
  sc_wrapper->bp_comparative        = NULL;
  sc_wrapper->bp_local              = NULL;
  sc_wrapper->bp_local_comparative  = NULL;
  sc_wrapper->stack                 = NULL;
  sc_wrapper->stack_comparative     = NULL;
  sc_wrapper->user_cb               = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data             = NULL;
  sc_wrapper->user_data_comparative = NULL;

  sc_wrapper->pair      = NULL;
  sc_wrapper->pair_ext  = NULL;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;

      if (sc) {
        sc_wrapper->up        = sc->exp_energy_up;
        sc_wrapper->bp        = (sliding_window) ? NULL : sc->exp_energy_bp;
        sc_wrapper->bp_local  = (sliding_window) ? sc->exp_energy_bp_local : NULL;
        sc_wrapper->stack     = sc->exp_energy_stack;
        sc_wrapper->user_cb   = sc->exp_f;
        sc_wrapper->user_data = sc->data;

        if (sc->exp_energy_up)
          provides_sc_up = 1;

        if (sliding_window) {
          if (sc->exp_energy_bp_local)
            provides_sc_bp = 1;
        } else if (sc->exp_energy_bp) {
          provides_sc_bp = 1;
        }

        if (sc->exp_energy_stack)
          provides_sc_stack = 1;

        if (sc->exp_f)
          provides_sc_user = 1;

        if (provides_sc_user) {
          if (provides_sc_up) {
            if (provides_sc_bp) {
              if (provides_sc_stack) {
                sc_wrapper->pair = (sliding_window) ?
                                   &sc_int_exp_cb_up_bp_local_stack_user :
                                   &sc_int_exp_cb_up_bp_stack_user;
                sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_stack_user;
              } else {
                sc_wrapper->pair = (sliding_window) ?
                                   &sc_int_exp_cb_up_bp_local_user :
                                   &sc_int_exp_cb_up_bp_user;
                sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_user;
              }
            } else if (provides_sc_stack) {
              sc_wrapper->pair      = &sc_int_exp_cb_up_stack_user;
              sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_stack_user;
            } else {
              sc_wrapper->pair      = &sc_int_exp_cb_up_user;
              sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_user;
            }
          } else if (provides_sc_bp) {
            if (provides_sc_stack) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_bp_local_stack_user :
                                 &sc_int_exp_cb_bp_stack_user;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_stack_user;
            } else {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_bp_local_user :
                                 &sc_int_exp_cb_bp_user;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_user;
            }
          } else if (provides_sc_stack) {
            sc_wrapper->pair      = &sc_int_exp_cb_stack_user;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_stack_user;
          } else {
            sc_wrapper->pair      = &sc_int_exp_cb_user;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_user;
          }
        } else if (provides_sc_bp) {
          if (provides_sc_up) {
            if (provides_sc_stack) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_up_bp_local_stack :
                                 &sc_int_exp_cb_up_bp_stack;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_stack;
            } else {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_up_bp_local :
                                 &sc_int_exp_cb_up_bp;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up;
            }
          } else if (provides_sc_stack) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_int_exp_cb_bp_local_stack :
                               &sc_int_exp_cb_bp_stack;
            sc_wrapper->pair_ext = &sc_int_exp_cb_ext_stack;
          } else {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_int_exp_cb_bp_local :
                               &sc_int_exp_cb_bp;
          }
        } else if (provides_sc_up) {
          if (provides_sc_stack) {
            sc_wrapper->pair      = &sc_int_exp_cb_up_stack;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_stack;
          } else {
            sc_wrapper->pair      = &sc_int_exp_cb_up;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up;
          }
        } else if (provides_sc_stack) {
          sc_wrapper->pair      = &sc_int_exp_cb_stack;
          sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_stack;
        }
      }

      break;

    case VRNA_FC_TYPE_COMPARATIVE:
      sc_wrapper->n_seq = fc->n_seq;
      sc_wrapper->a2s   = fc->a2s;
      scs               = fc->scs;

      if (scs) {
        sc_wrapper->up_comparative = (FLT_OR_DBL ***)vrna_alloc(
          sizeof(FLT_OR_DBL * *) *
          fc->n_seq);

        sc_wrapper->bp_comparative = (FLT_OR_DBL **)vrna_alloc(
          sizeof(FLT_OR_DBL *) *
          fc->n_seq);

        sc_wrapper->bp_local_comparative = (FLT_OR_DBL ***)vrna_alloc(
          sizeof(FLT_OR_DBL * *) *
          fc->n_seq);

        sc_wrapper->stack_comparative = (FLT_OR_DBL **)vrna_alloc(
          sizeof(FLT_OR_DBL *) *
          fc->n_seq);

        sc_wrapper->user_cb_comparative = (vrna_sc_exp_f *)vrna_alloc(
          sizeof(vrna_sc_exp_f) *
          fc->n_seq);

        sc_wrapper->user_data_comparative = (void **)vrna_alloc(
          sizeof(void *) *
          fc->n_seq);

        for (s = 0; s < fc->n_seq; s++) {
          if (scs[s]) {
            sliding_window = (scs[s]->type == VRNA_SC_WINDOW) ?
                             1 :
                             0;
            sc_wrapper->up_comparative[s] = scs[s]->exp_energy_up;
            sc_wrapper->bp_comparative[s] = (sliding_window) ?
                                            NULL :
                                            scs[s]->exp_energy_bp;
            sc_wrapper->bp_local_comparative[s] = (sliding_window) ?
                                                  scs[s]->exp_energy_bp_local :
                                                  NULL;
            sc_wrapper->stack_comparative[s]      = scs[s]->exp_energy_stack;
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

            if (scs[s]->exp_energy_stack)
              provides_sc_stack = 1;

            if (scs[s]->exp_f)
              provides_sc_user = 1;
          }
        }

        if (provides_sc_user) {
          if (provides_sc_up) {
            if (provides_sc_bp) {
              if (provides_sc_stack) {
                sc_wrapper->pair = (sliding_window) ?
                                   &sc_int_exp_cb_up_bp_local_stack_user_comparative :
                                   &sc_int_exp_cb_up_bp_stack_user_comparative;
                sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_stack_user_comparative;
              } else {
                sc_wrapper->pair = (sliding_window) ?
                                   &sc_int_exp_cb_up_bp_local_user_comparative :
                                   &sc_int_exp_cb_up_bp_user_comparative;
                sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_user_comparative;
              }
            } else if (provides_sc_stack) {
              sc_wrapper->pair      = &sc_int_exp_cb_up_stack_user_comparative;
              sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_stack_user_comparative;
            } else {
              sc_wrapper->pair      = &sc_int_exp_cb_up_user_comparative;
              sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_user_comparative;
            }
          } else if (provides_sc_bp) {
            if (provides_sc_stack) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_bp_local_stack_user_comparative :
                                 &sc_int_exp_cb_bp_stack_user_comparative;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_stack_user_comparative;
            } else {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_bp_local_user_comparative :
                                 &sc_int_exp_cb_bp_user_comparative;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_user_comparative;
            }
          } else if (provides_sc_stack) {
            sc_wrapper->pair      = &sc_int_exp_cb_stack_user_comparative;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_stack_user_comparative;
          } else {
            sc_wrapper->pair      = &sc_int_exp_cb_user_comparative;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_user_comparative;
          }
        } else if (provides_sc_bp) {
          if (provides_sc_up) {
            if (provides_sc_stack) {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_up_bp_local_stack_comparative :
                                 &sc_int_exp_cb_up_bp_stack_comparative;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_stack_comparative;
            } else {
              sc_wrapper->pair = (sliding_window) ?
                                 &sc_int_exp_cb_up_bp_local_comparative :
                                 &sc_int_exp_cb_up_bp_comparative;
              sc_wrapper->pair_ext = &sc_int_exp_cb_ext_up_comparative;
            }
          } else if (provides_sc_stack) {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_int_exp_cb_bp_local_stack_comparative :
                               &sc_int_exp_cb_bp_stack_comparative;
            sc_wrapper->pair_ext = &sc_int_exp_cb_ext_stack_comparative;
          } else {
            sc_wrapper->pair = (sliding_window) ?
                               &sc_int_exp_cb_bp_local_comparative :
                               &sc_int_exp_cb_bp_comparative;
          }
        } else if (provides_sc_up) {
          if (provides_sc_stack) {
            sc_wrapper->pair      = &sc_int_exp_cb_up_stack_comparative;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_stack_comparative;
          } else {
            sc_wrapper->pair      = &sc_int_exp_cb_up_comparative;
            sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_up_comparative;
          }
        } else if (provides_sc_stack) {
          sc_wrapper->pair      = &sc_int_exp_cb_stack_comparative;
          sc_wrapper->pair_ext  = &sc_int_exp_cb_ext_stack_comparative;
        }
      }

      break;
  }
}


PRIVATE INLINE void
free_sc_int_exp(struct sc_int_exp_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->bp_comparative);
  free(sc_wrapper->bp_local_comparative);
  free(sc_wrapper->stack_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}
