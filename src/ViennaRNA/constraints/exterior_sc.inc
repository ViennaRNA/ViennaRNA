/*
 *  This file contains soft constraints wrappers for exterior loop evaluation
 */

#include <ViennaRNA/utils/basic.h>

struct sc_f5_dat;
struct sc_f3_dat;

typedef int (*sc_f5_cb)(int              j,
                       int              k,
                       int              l,
                       struct sc_f5_dat *data);


typedef int (*sc_f3_cb)(int              j,
                       int              k,
                       int              l,
                       struct sc_f3_dat *data);


typedef int (*sc_ext_pair_cb)(int               i,
                             int               j,
                             struct sc_f5_dat  *data);

typedef int (*sc_ext_red_cb)(int              i,
                            int              j,
                            int              k,
                            int              l,
                            struct sc_f5_dat *data);

struct sc_f5_dat {
  int                     **up;

  sc_f5_cb                red_ext5;
  sc_f5_cb                red_stem5;
  sc_f5_cb                decomp_stem5;
  sc_f5_cb                decomp_stem51;

  sc_ext_pair_cb          pair;
  sc_ext_red_cb           red_ext;
  sc_ext_red_cb           red_stem;
  sc_ext_red_cb           decomp;

  vrna_sc_f user_cb;
  void                    *user_data;

  /* below attributes are for comparative structure prediction */
  unsigned int            n_seq;
  unsigned int            **a2s;
  int                     ***up_comparative;

  vrna_sc_f *user_cb_comparative;
  void                    **user_data_comparative;
};


struct sc_f3_dat {
  unsigned int            n;

  int                     **up;

  sc_f3_cb                red_ext;
  sc_f3_cb                red_stem;
  sc_f3_cb                decomp_stem;
  sc_f3_cb                decomp_stem1;

  vrna_sc_f user_cb;
  void                    *user_data;

  /* below attributes are for comparative structure prediction */
  unsigned int            n_seq;
  unsigned int            **a2s;
  int                     ***up_comparative;

  vrna_sc_f *user_cb_comparative;
  void                    **user_data_comparative;
};


PRIVATE INLINE void
init_sc_f5(vrna_fold_compound_t *fc,
           struct sc_f5_dat     *sc_wrapper);


PRIVATE INLINE void
init_sc_f3(vrna_fold_compound_t *fc,
           int                  i,
           struct sc_f3_dat     *sc_wrapper);


PRIVATE INLINE void
free_sc_f5(struct sc_f5_dat *sc_wrapper);


PRIVATE INLINE void
free_sc_f3(struct sc_f3_dat *sc_wrapper);


PRIVATE INLINE int
sc_f5_cb_reduce(int               j,
                int               k,
                int               l,
                struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_ext_cb_reduce(int i,
                int               j,
                int               k,
                int               l,
                struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_reduce(int               i,
                int               k,
                int               l,
                struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_reduce_comparative(int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_ext_cb_reduce_comparative(int i,
                            int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_reduce_comparative(int               i,
                            int               k,
                            int               l,
                            struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_def_reduce_to_ext(int               j,
                                int               k,
                                int               l,
                                struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_ext_cb_user_def_reduce_to_ext(int i,
                                int               j,
                                int               k,
                                int               l,
                                struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_user_def_reduce_to_ext(int               i,
                                int               k,
                                int               l,
                                struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_reduce_to_ext(int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data);

PRIVATE INLINE int
sc_ext_cb_user_reduce_to_ext(int i,
                              int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_user_reduce_to_ext(int               i,
                            int               k,
                            int               l,
                            struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_def_reduce_to_ext_comparative(int               j,
                                            int               k,
                                            int               l,
                                            struct sc_f5_dat  *data);

PRIVATE INLINE int
sc_ext_cb_user_def_reduce_to_ext_comparative(int i,
                                              int               j,
                                            int               k,
                                            int               l,
                                            struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_user_def_reduce_to_ext_comparative(int               i,
                                            int               k,
                                            int               l,
                                            struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_reduce_to_ext_comparative(int               j,
                                        int               k,
                                        int               l,
                                        struct sc_f5_dat  *data);

PRIVATE INLINE int
sc_ext_cb_user_reduce_to_ext_comparative(int i,
                                        int               j,
                                        int               k,
                                        int               l,
                                        struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_user_reduce_to_ext_comparative(int               i,
                                        int               k,
                                        int               l,
                                        struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_def_reduce_to_stem(int              j,
                                 int              k,
                                 int              l,
                                 struct sc_f5_dat *data);

PRIVATE INLINE int
sc_ext_cb_user_def_reduce_to_stem(int i,
                                  int              j,
                                 int              k,
                                 int              l,
                                 struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_user_def_reduce_to_stem(int              i,
                                 int              k,
                                 int              l,
                                 struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_user_reduce_to_stem(int              j,
                             int              k,
                             int              l,
                             struct sc_f5_dat *data);

PRIVATE INLINE int
sc_ext_cb_user_reduce_to_stem(int i,
                                int              j,
                             int              k,
                             int              l,
                             struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_user_reduce_to_stem(int              i,
                             int              k,
                             int              l,
                             struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_user_def_reduce_to_stem_comparative(int              j,
                                             int              k,
                                             int              l,
                                             struct sc_f5_dat *data);

PRIVATE INLINE int
sc_ext_cb_user_def_reduce_to_stem_comparative(int i,
                                            int              j,
                                             int              k,
                                             int              l,
                                             struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_user_def_reduce_to_stem_comparative(int              i,
                                             int              k,
                                             int              l,
                                             struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_user_reduce_to_stem_comparative(int              j,
                                         int              k,
                                         int              l,
                                         struct sc_f5_dat *data);

PRIVATE INLINE int
sc_ext_cb_user_reduce_to_stem_comparative(int i,
                                          int              j,
                                         int              k,
                                         int              l,
                                         struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_user_reduce_to_stem_comparative(int              i,
                                         int              k,
                                         int              l,
                                         struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_split_in_ext_stem(int              j,
                           int              k,
                           int              l,
                           struct sc_f5_dat *data);

PRIVATE INLINE int
sc_ext_cb_split(int i,
                            int              j,
                           int              k,
                           int              l,
                           struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_split_in_stem_ext(int              i,
                           int              k,
                           int              l,
                           struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_user_def_split_in_ext_stem(int               j,
                                    int               k,
                                    int               l,
                                    struct sc_f5_dat  *data);

PRIVATE INLINE int
sc_ext_cb_user_def_split(int i,
                                      int               j,
                                    int               k,
                                    int               l,
                                    struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_user_def_split_in_stem_ext(int               i,
                                    int               k,
                                    int               l,
                                    struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_split_in_ext_stem(int               j,
                                int               k,
                                int               l,
                                struct sc_f5_dat  *data);

PRIVATE INLINE int
sc_ext_cb_user_split(int i,
                                int               j,
                                int               k,
                                int               l,
                                struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_user_split_in_stem_ext(int               i,
                                int               k,
                                int               l,
                                struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_split_in_ext_stem_comparative(int              j,
                                       int              k,
                                       int              l,
                                       struct sc_f5_dat *data);

PRIVATE INLINE int
sc_ext_cb_split_comparative(int i,
                                        int              j,
                                       int              k,
                                       int              l,
                                       struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_split_in_stem_ext_comparative(int              i,
                                       int              k,
                                       int              l,
                                       struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_user_def_split_in_ext_stem_comparative(int               j,
                                                int               k,
                                                int               l,
                                                struct sc_f5_dat  *data);

PRIVATE INLINE int
sc_ext_cb_user_def_split_comparative(int i,
                                                int               j,
                                                int               k,
                                                int               l,
                                                struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_user_def_split_in_stem_ext_comparative(int               i,
                                                int               k,
                                                int               l,
                                                struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_split_in_ext_stem_comparative(int               j,
                                            int               k,
                                            int               l,
                                            struct sc_f5_dat  *data);

PRIVATE INLINE int
sc_ext_cb_user_split_comparative(int i,
                                            int               j,
                                            int               k,
                                            int               l,
                                            struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_user_split_in_stem_ext_comparative(int               i,
                                            int               k,
                                            int               l,
                                            struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_split_in_ext_stem1(int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_split_in_stem_ext1(int               i,
                            int               k,
                            int               l,
                            struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_def_split_in_ext_stem1(int              j,
                                     int              k,
                                     int              l,
                                     struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_user_def_split_in_stem_ext1(int              i,
                                     int              k,
                                     int              l,
                                     struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_user_split_in_ext_stem1(int              j,
                                 int              k,
                                 int              l,
                                 struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_user_split_in_stem_ext1(int              i,
                                 int              k,
                                 int              l,
                                 struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_split_in_ext_stem1_comparative(int               j,
                                        int               k,
                                        int               l,
                                        struct sc_f5_dat  *data);


PRIVATE INLINE int
sc_f3_cb_split_in_stem_ext1_comparative(int               i,
                                        int               k,
                                        int               l,
                                        struct sc_f3_dat  *data);


PRIVATE INLINE int
sc_f5_cb_user_def_split_in_ext_stem1_comparative(int              j,
                                                 int              k,
                                                 int              l,
                                                 struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_user_def_split_in_stem_ext1_comparative(int              i,
                                                 int              k,
                                                 int              l,
                                                 struct sc_f3_dat *data);


PRIVATE INLINE int
sc_f5_cb_user_split_in_ext_stem1_comparative(int              j,
                                             int              k,
                                             int              l,
                                             struct sc_f5_dat *data);


PRIVATE INLINE int
sc_f3_cb_user_split_in_stem_ext1_comparative(int              i,
                                             int              k,
                                             int              l,
                                             struct sc_f3_dat *data);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PRIVATE INLINE void
init_sc_f5(vrna_fold_compound_t *fc,
           struct sc_f5_dat     *sc_wrapper)
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
  sc_wrapper->red_ext5      = NULL;
  sc_wrapper->red_stem5     = NULL;
  sc_wrapper->decomp_stem5  = NULL;
  sc_wrapper->decomp_stem51 = NULL;

  sc_wrapper->pair= NULL;
  sc_wrapper->red_ext= NULL;
  sc_wrapper->red_stem= NULL;
  sc_wrapper->decomp = NULL;

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
            sc_wrapper->red_ext5      = &sc_f5_cb_user_def_reduce_to_ext;
            sc_wrapper->red_stem5      = &sc_f5_cb_user_def_reduce_to_stem;
            sc_wrapper->decomp_stem5   = &sc_f5_cb_user_def_split_in_ext_stem;
            sc_wrapper->decomp_stem51  = &sc_f5_cb_user_def_split_in_ext_stem1;
            sc_wrapper->red_ext         = &sc_ext_cb_user_def_reduce_to_ext;
            sc_wrapper->red_stem      = &sc_ext_cb_user_def_reduce_to_stem;
            sc_wrapper->decomp   = &sc_ext_cb_user_def_split;
          } else {
            sc_wrapper->red_ext5       = &sc_f5_cb_reduce;
            sc_wrapper->red_stem5      = &sc_f5_cb_reduce;
            sc_wrapper->decomp_stem5   = &sc_f5_cb_split_in_ext_stem;
            sc_wrapper->decomp_stem51  = &sc_f5_cb_split_in_ext_stem1;
            sc_wrapper->red_ext       = &sc_ext_cb_reduce;
            sc_wrapper->red_stem      = &sc_ext_cb_reduce;
            sc_wrapper->decomp   = &sc_ext_cb_split;
          }
        } else if (sc->f) {
          sc_wrapper->red_ext5       = &sc_f5_cb_user_reduce_to_ext;
          sc_wrapper->red_stem5      = &sc_f5_cb_user_reduce_to_stem;
          sc_wrapper->decomp_stem5   = &sc_f5_cb_user_split_in_ext_stem;
          sc_wrapper->decomp_stem51  = &sc_f5_cb_user_split_in_ext_stem1;
          sc_wrapper->red_ext       = &sc_ext_cb_user_reduce_to_ext;
          sc_wrapper->red_stem      = &sc_ext_cb_user_reduce_to_stem;
          sc_wrapper->decomp   = &sc_ext_cb_user_split;
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
          (vrna_sc_f *)vrna_alloc(sizeof(vrna_sc_f) * fc->n_seq);
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
            sc_wrapper->red_ext5       = &sc_f5_cb_user_def_reduce_to_ext_comparative;
            sc_wrapper->red_stem5      = &sc_f5_cb_user_def_reduce_to_stem_comparative;
            sc_wrapper->decomp_stem5   = &sc_f5_cb_user_def_split_in_ext_stem_comparative;
            sc_wrapper->decomp_stem51  = &sc_f5_cb_user_def_split_in_ext_stem1_comparative;
            sc_wrapper->red_ext       = &sc_ext_cb_user_def_reduce_to_ext_comparative;
            sc_wrapper->red_stem      = &sc_ext_cb_user_def_reduce_to_stem_comparative;
            sc_wrapper->decomp   = &sc_ext_cb_user_def_split;
          } else {
            sc_wrapper->red_ext5       = &sc_f5_cb_reduce_comparative;
            sc_wrapper->red_stem5      = &sc_f5_cb_reduce_comparative;
            sc_wrapper->decomp_stem5   = &sc_f5_cb_split_in_ext_stem_comparative;
            sc_wrapper->decomp_stem51  = &sc_f5_cb_split_in_ext_stem1_comparative;
            sc_wrapper->red_ext       = &sc_ext_cb_reduce_comparative;
            sc_wrapper->red_stem      = &sc_ext_cb_reduce_comparative;
            sc_wrapper->decomp   = &sc_ext_cb_split;
          }
        } else if (provides_sc_user_cb) {
          sc_wrapper->red_ext5       = &sc_f5_cb_user_reduce_to_ext_comparative;
          sc_wrapper->red_stem5      = &sc_f5_cb_user_reduce_to_stem_comparative;
          sc_wrapper->decomp_stem5   = &sc_f5_cb_user_split_in_ext_stem_comparative;
          sc_wrapper->decomp_stem51  = &sc_f5_cb_user_split_in_ext_stem1_comparative;
          sc_wrapper->red_ext       = &sc_ext_cb_user_reduce_to_ext_comparative;
          sc_wrapper->red_stem      = &sc_ext_cb_user_reduce_to_stem_comparative;
          sc_wrapper->decomp   = &sc_ext_cb_user_split;
        }
      }

      break;
  }
}


PRIVATE INLINE void
free_sc_f5(struct sc_f5_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}


static void
init_sc_f3(vrna_fold_compound_t *fc,
           int                  i VRNA_UNUSED,
           struct sc_f3_dat     *sc_wrapper)
{
  unsigned int  s;
  vrna_sc_t     *sc, **scs;

  sc_wrapper->n             = fc->length;
  sc_wrapper->n_seq         = 1;
  sc_wrapper->up            = NULL;
  sc_wrapper->red_ext       = NULL;
  sc_wrapper->red_stem      = NULL;
  sc_wrapper->decomp_stem   = NULL;
  sc_wrapper->decomp_stem1  = NULL;
  sc_wrapper->user_cb       = NULL;
  sc_wrapper->user_data     = NULL;

  sc_wrapper->up_comparative        = NULL;
  sc_wrapper->user_cb_comparative   = NULL;
  sc_wrapper->user_data_comparative = NULL;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sc = fc->sc;
      if (sc) {
        sc_wrapper->up        = sc->energy_up;
        sc_wrapper->user_cb   = sc->f;
        sc_wrapper->user_data = sc->data;

        if (sc->energy_up) {
          if (sc->f) {
            sc_wrapper->red_ext       = &sc_f3_cb_user_def_reduce_to_ext;
            sc_wrapper->red_stem      = &sc_f3_cb_user_def_reduce_to_stem;
            sc_wrapper->decomp_stem   = &sc_f3_cb_user_def_split_in_stem_ext;
            sc_wrapper->decomp_stem1  = &sc_f3_cb_user_def_split_in_stem_ext1;
          } else {
            sc_wrapper->red_ext       = &sc_f3_cb_reduce;
            sc_wrapper->red_stem      = &sc_f3_cb_reduce;
            sc_wrapper->decomp_stem   = &sc_f3_cb_split_in_stem_ext;
            sc_wrapper->decomp_stem1  = &sc_f3_cb_split_in_stem_ext1;
          }
        } else if (sc->f) {
          sc_wrapper->red_ext       = &sc_f3_cb_user_reduce_to_ext;
          sc_wrapper->red_stem      = &sc_f3_cb_user_reduce_to_stem;
          sc_wrapper->decomp_stem   = &sc_f3_cb_user_split_in_stem_ext;
          sc_wrapper->decomp_stem1  = &sc_f3_cb_user_split_in_stem_ext1;
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
          (vrna_sc_f *)vrna_alloc(sizeof(vrna_sc_f) * fc->n_seq);
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
            sc_wrapper->red_ext       = &sc_f3_cb_user_def_reduce_to_ext_comparative;
            sc_wrapper->red_stem      = &sc_f3_cb_user_def_reduce_to_stem_comparative;
            sc_wrapper->decomp_stem   = &sc_f3_cb_user_def_split_in_stem_ext_comparative;
            sc_wrapper->decomp_stem1  = &sc_f3_cb_user_def_split_in_stem_ext1_comparative;
          } else {
            sc_wrapper->red_ext       = &sc_f3_cb_reduce_comparative;
            sc_wrapper->red_stem      = &sc_f3_cb_reduce_comparative;
            sc_wrapper->decomp_stem   = &sc_f3_cb_split_in_stem_ext_comparative;
            sc_wrapper->decomp_stem1  = &sc_f3_cb_split_in_stem_ext1_comparative;
          }
        } else if (provides_sc_user_cb) {
          sc_wrapper->red_ext       = &sc_f3_cb_user_reduce_to_ext_comparative;
          sc_wrapper->red_stem      = &sc_f3_cb_user_reduce_to_stem_comparative;
          sc_wrapper->decomp_stem   = &sc_f3_cb_user_split_in_stem_ext_comparative;
          sc_wrapper->decomp_stem1  = &sc_f3_cb_user_split_in_stem_ext1_comparative;
        }
      }

      break;
  }
}


static void
free_sc_f3(struct sc_f3_dat *sc_wrapper)
{
  free(sc_wrapper->up_comparative);
  free(sc_wrapper->user_cb_comparative);
  free(sc_wrapper->user_data_comparative);
}


PRIVATE INLINE int
sc_f5_cb_reduce(int               j,
                int               k,
                int               l,
                struct sc_f5_dat  *data)
{
  return sc_ext_cb_reduce(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_reduce(int i,
                int               j,
                int               k,
                int               l,
                struct sc_f5_dat  *data)
{
  unsigned int  start_2, length_1, length_2;
  int           e_sc, **sc_up;

  sc_up = data->up;

  e_sc = 0;

  length_1  = k - i;
  start_2   = l + 1;
  length_2  = j - l;

  if (length_1 != 0)
    e_sc += sc_up[i][length_1];

  if (length_2 != 0)
    e_sc += sc_up[start_2][length_2];

  return e_sc;
}


PRIVATE INLINE int
sc_f3_cb_reduce(int               i,
                int               k,
                int               l,
                struct sc_f3_dat  *data)
{
  unsigned int  start_2, length_1, length_2;
  int           e_sc, **sc_up;

  sc_up = data->up;

  e_sc = 0;

  length_1  = k - i;
  start_2   = l + 1;
  length_2  = data->n - l;

  if (length_1 != 0)
    e_sc += sc_up[i][length_1];

  if (length_2 != 0)
    e_sc += sc_up[start_2][length_2];

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_reduce_comparative(int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data)
{
  return sc_ext_cb_reduce_comparative(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_reduce_comparative(int              i,
                             int              j,
                             int              k,
                             int              l,
                             struct sc_f5_dat *data)
{
  unsigned int  s, start_1, start_2, length_1, length_2, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;
  a2s   = data->a2s;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][i];
      length_1  = a2s[s][k] - start_1;
      start_2   = a2s[s][l] + 1;
      length_2  = a2s[s][j] - a2s[s][l];

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];

      if (length_2 != 0)
        e_sc += sc_up[s][start_2][length_2];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_f3_cb_reduce_comparative(int               i,
                            int               k,
                            int               l,
                            struct sc_f3_dat  *data)
{
  unsigned int  s, start_1, start_2, length_1, length_2, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;
  a2s   = data->a2s;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][i];
      length_1  = a2s[s][k] - start_1;
      start_2   = a2s[s][l] + 1;
      length_2  = a2s[s][data->n] - a2s[s][l];

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];

      if (length_2 != 0)
        e_sc += sc_up[s][start_2][length_2];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_user_reduce_to_ext(int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data)
{
  return data->user_cb(1, j, k, l, VRNA_DECOMP_EXT_EXT, data->user_data);
}


PRIVATE INLINE int
sc_ext_cb_user_reduce_to_ext(int i,
                            int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data)
{
  return data->user_cb(i, j, k, l, VRNA_DECOMP_EXT_EXT, data->user_data);
}


PRIVATE INLINE int
sc_f3_cb_user_reduce_to_ext(int               i,
                            int               k,
                            int               l,
                            struct sc_f3_dat  *data)
{
  return data->user_cb(i, data->n, k, l, VRNA_DECOMP_EXT_EXT, data->user_data);
}


PRIVATE INLINE int
sc_f5_cb_user_reduce_to_ext_comparative(int               j,
                                        int               k,
                                        int               l,
                                        struct sc_f5_dat  *data)
{
  return sc_ext_cb_user_reduce_to_ext_comparative(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_user_reduce_to_ext_comparative(int i,
                                          int               j,
                                        int               k,
                                        int               l,
                                        struct sc_f5_dat  *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc += data->user_cb_comparative[s](i,
                                           j,
                                           k,
                                           l,
                                           VRNA_DECOMP_EXT_EXT,
                                           data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_f3_cb_user_reduce_to_ext_comparative(int               i,
                                        int               k,
                                        int               l,
                                        struct sc_f3_dat  *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc +=
        data->user_cb_comparative[s](i, data->n, k, l, VRNA_DECOMP_EXT_EXT,
                                     data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_user_def_reduce_to_ext(int               j,
                                int               k,
                                int               l,
                                struct sc_f5_dat  *data)
{
  return sc_ext_cb_user_def_reduce_to_ext(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_user_def_reduce_to_ext(int i,
                                int               j,
                                int               k,
                                int               l,
                                struct sc_f5_dat  *data)
{
  return sc_ext_cb_reduce(i, j, k, l, data) +
         sc_ext_cb_user_reduce_to_ext(i, j, k, l, data);
}


PRIVATE INLINE int
sc_f3_cb_user_def_reduce_to_ext(int               i,
                                int               k,
                                int               l,
                                struct sc_f3_dat  *data)
{
  return sc_f3_cb_reduce(i, k, l, data) +
         sc_f3_cb_user_reduce_to_ext(i, k, l, data);
}


PRIVATE INLINE int
sc_f5_cb_user_def_reduce_to_ext_comparative(int               j,
                                            int               k,
                                            int               l,
                                            struct sc_f5_dat  *data)
{
  return sc_ext_cb_user_def_reduce_to_ext_comparative(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_user_def_reduce_to_ext_comparative(int i,
                                              int               j,
                                            int               k,
                                            int               l,
                                            struct sc_f5_dat  *data)
{
  return sc_ext_cb_reduce_comparative(i, j, k, l, data) +
         sc_ext_cb_user_reduce_to_ext_comparative(i, j, k, l, data);
}


PRIVATE INLINE int
sc_f3_cb_user_def_reduce_to_ext_comparative(int               i,
                                            int               k,
                                            int               l,
                                            struct sc_f3_dat  *data)
{
  return sc_f3_cb_reduce_comparative(i, k, l, data) +
         sc_f3_cb_user_reduce_to_ext_comparative(i, k, l, data);
}


PRIVATE INLINE int
sc_f5_cb_user_reduce_to_stem(int              j,
                             int              k,
                             int              l,
                             struct sc_f5_dat *data)
{
  return sc_ext_cb_user_reduce_to_stem(1, j, k, l, data);
}

PRIVATE INLINE int
sc_ext_cb_user_reduce_to_stem(int i,
                            int              j,
                             int              k,
                             int              l,
                             struct sc_f5_dat *data)
{
  return data->user_cb(i, j, k, l, VRNA_DECOMP_EXT_STEM, data->user_data);
}


PRIVATE INLINE int
sc_f3_cb_user_reduce_to_stem(int              i,
                             int              k,
                             int              l,
                             struct sc_f3_dat *data)
{
  return data->user_cb(i, data->n, k, l, VRNA_DECOMP_EXT_STEM, data->user_data);
}


PRIVATE INLINE int
sc_f5_cb_user_reduce_to_stem_comparative(int              j,
                                         int              k,
                                         int              l,
                                         struct sc_f5_dat *data)
{
  return sc_ext_cb_user_reduce_to_stem_comparative(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_user_reduce_to_stem_comparative(int i,
                                          int              j,
                                         int              k,
                                         int              l,
                                         struct sc_f5_dat *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc += data->user_cb_comparative[s](i,
                                           j,
                                           k,
                                           l,
                                           VRNA_DECOMP_EXT_STEM,
                                           data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_f3_cb_user_reduce_to_stem_comparative(int              i,
                                         int              k,
                                         int              l,
                                         struct sc_f3_dat *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc +=
        data->user_cb_comparative[s](i, data->n, k, l, VRNA_DECOMP_EXT_STEM,
                                     data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_user_def_reduce_to_stem(int              j,
                                 int              k,
                                 int              l,
                                 struct sc_f5_dat *data)
{
  return sc_ext_cb_user_def_reduce_to_stem(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_user_def_reduce_to_stem(int i,
                                  int              j,
                                 int              k,
                                 int              l,
                                 struct sc_f5_dat *data)
{
  return sc_ext_cb_reduce(i, j, k, l, data) +
         sc_ext_cb_user_reduce_to_stem(i, j, k, l, data);
}


PRIVATE INLINE int
sc_f3_cb_user_def_reduce_to_stem(int              i,
                                 int              k,
                                 int              l,
                                 struct sc_f3_dat *data)
{
  return sc_f3_cb_reduce(i, k, l, data) +
         sc_f3_cb_user_reduce_to_stem(i, k, l, data);
}


PRIVATE INLINE int
sc_f5_cb_user_def_reduce_to_stem_comparative(int              j,
                                             int              k,
                                             int              l,
                                             struct sc_f5_dat *data)
{
  return sc_ext_cb_user_def_reduce_to_stem_comparative(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_user_def_reduce_to_stem_comparative(int i,
                                              int              j,
                                             int              k,
                                             int              l,
                                             struct sc_f5_dat *data)
{
  return sc_ext_cb_reduce_comparative(i, j, k, l, data) +
         sc_ext_cb_user_reduce_to_stem_comparative(i, j, k, l, data);
}


PRIVATE INLINE int
sc_f3_cb_user_def_reduce_to_stem_comparative(int              i,
                                             int              k,
                                             int              l,
                                             struct sc_f3_dat *data)
{
  return sc_f3_cb_reduce_comparative(i, k, l, data) +
         sc_f3_cb_user_reduce_to_stem_comparative(i, k, l, data);
}


PRIVATE INLINE int
sc_f5_cb_split_in_ext_stem(int              j,
                           int              k,
                           int              l,
                           struct sc_f5_dat *data)
{
  return sc_ext_cb_split(1, j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_split(int              i VRNA_UNUSED,
                int              j VRNA_UNUSED,
                int              k,
                int              l,
                struct sc_f5_dat *data)
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
sc_f3_cb_split_in_stem_ext(int              i VRNA_UNUSED,
                           int              k,
                           int              l,
                           struct sc_f3_dat *data)
{
  unsigned int  start_3, length_3;
  int           e_sc, **sc_up;

  sc_up = data->up;

  e_sc = 0;

  start_3   = k + 1;
  length_3  = l - k - 1;

  if (length_3 != 0)
    e_sc += sc_up[start_3][length_3];

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_split_in_ext_stem_comparative(int              j VRNA_UNUSED,
                                       int              k,
                                       int              l,
                                       struct sc_f5_dat *data)
{
  unsigned int  s, start_1, length_1, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;
  a2s   = data->a2s;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][k] + 1;
      length_1  = a2s[s][l - 1] - a2s[s][k];

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];
    }
  }

  return e_sc;
}

PRIVATE INLINE int
sc_ext_cb_split_comparative(int               i VRNA_UNUSED,
                            int               j VRNA_UNUSED,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data)
{
  unsigned int  s, start_1, length_1, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;
  a2s   = data->a2s;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][k] + 1;
      length_1  = a2s[s][l - 1] - a2s[s][k];

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_f3_cb_split_in_stem_ext_comparative(int              i VRNA_UNUSED,
                                       int              k,
                                       int              l,
                                       struct sc_f3_dat *data)
{
  unsigned int  s, start_1, length_1, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;
  a2s   = data->a2s;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][k] + 1;
      length_1  = a2s[s][l - 1] - a2s[s][k];

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_user_split_in_ext_stem(int               j,
                                int               k,
                                int               l,
                                struct sc_f5_dat  *data)
{
  return data->user_cb(1, j, k, l, VRNA_DECOMP_EXT_EXT_STEM, data->user_data);
}


PRIVATE INLINE int
sc_ext_cb_user_split(int i,
                                int               j,
                                int               k,
                                int               l,
                                struct sc_f5_dat  *data)
{
  return data->user_cb(i, j, k, l, VRNA_DECOMP_EXT_EXT_EXT, data->user_data);
}


PRIVATE INLINE int
sc_f3_cb_user_split_in_stem_ext(int               i,
                                int               k,
                                int               l,
                                struct sc_f3_dat  *data)
{
  return data->user_cb(i, data->n, k, l, VRNA_DECOMP_EXT_STEM_EXT, data->user_data);
}


PRIVATE INLINE int
sc_f5_cb_user_split_in_ext_stem_comparative(int               j,
                                            int               k,
                                            int               l,
                                            struct sc_f5_dat  *data)
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
sc_ext_cb_user_split_comparative(int i,
                                            int               j,
                                            int               k,
                                            int               l,
                                            struct sc_f5_dat  *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc += data->user_cb_comparative[s](i,
                                           j,
                                           k,
                                           l,
                                           VRNA_DECOMP_EXT_EXT_EXT,
                                           data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_f3_cb_user_split_in_stem_ext_comparative(int               i,
                                            int               k,
                                            int               l,
                                            struct sc_f3_dat  *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc +=
        data->user_cb_comparative[s](i, data->n, k, l, VRNA_DECOMP_EXT_STEM_EXT,
                                     data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_user_def_split_in_ext_stem(int               j,
                                    int               k,
                                    int               l,
                                    struct sc_f5_dat  *data)
{
  return sc_f5_cb_split_in_ext_stem(j, k, l, data) +
         sc_f5_cb_user_split_in_ext_stem(j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_user_def_split(int i,
                                    int               j,
                                    int               k,
                                    int               l,
                                    struct sc_f5_dat  *data)
{
  return sc_ext_cb_split(i, j, k, l, data) +
         sc_ext_cb_user_split(i, j, k, l, data);
}


PRIVATE INLINE int
sc_f3_cb_user_def_split_in_stem_ext(int               i,
                                    int               k,
                                    int               l,
                                    struct sc_f3_dat  *data)
{
  return sc_f3_cb_split_in_stem_ext(i, k, l, data) +
         sc_f3_cb_user_split_in_stem_ext(i, k, l, data);
}


PRIVATE INLINE int
sc_f5_cb_user_def_split_in_ext_stem_comparative(int               j,
                                                int               k,
                                                int               l,
                                                struct sc_f5_dat  *data)
{
  return sc_f5_cb_split_in_ext_stem_comparative(j, k, l, data) +
         sc_f5_cb_user_split_in_ext_stem_comparative(j, k, l, data);
}


PRIVATE INLINE int
sc_ext_cb_user_def_split_comparative(int i,
                                                  int               j,
                                                int               k,
                                                int               l,
                                                struct sc_f5_dat  *data)
{
  return sc_ext_cb_split_comparative(i, j, k, l, data) +
         sc_ext_cb_user_split_comparative(i, j, k, l, data);
}


PRIVATE INLINE int
sc_f3_cb_user_def_split_in_stem_ext_comparative(int               i,
                                                int               k,
                                                int               l,
                                                struct sc_f3_dat  *data)
{
  return sc_f3_cb_split_in_stem_ext_comparative(i, k, l, data) +
         sc_f3_cb_user_split_in_stem_ext_comparative(i, k, l, data);
}


PRIVATE INLINE int
sc_f5_cb_split_in_ext_stem1(int               j,
                            int               k,
                            int               l,
                            struct sc_f5_dat  *data)
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
sc_f3_cb_split_in_stem_ext1(int               i,
                            int               k,
                            int               l,
                            struct sc_f3_dat  *data)
{
  unsigned int  start_1, start_2, length_1, length_2;
  int           e_sc, **sc_up;

  sc_up = data->up;

  e_sc = 0;

  start_1   = k + 1;
  length_1  = l - k - 1;
  start_2   = i;
  length_2  = 1;

  if (length_1 != 0)
    e_sc += sc_up[start_1][length_1];

  e_sc += sc_up[start_2][length_2];

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_split_in_ext_stem1_comparative(int               j,
                                        int               k,
                                        int               l,
                                        struct sc_f5_dat  *data)
{
  unsigned int  s, start_1, start_2, length_1, length_2, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;

  e_sc  = 0;
  a2s   = data->a2s;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][k] + 1;
      length_1  = a2s[s][l - 1] - a2s[s][k];
      start_2   = a2s[s][j - 1];
      length_2  = a2s[s][j] - a2s[s][j - 1];

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];

      if (length_2 != 0)
        e_sc += sc_up[s][start_2][length_2];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_f3_cb_split_in_stem_ext1_comparative(int               i,
                                        int               k,
                                        int               l,
                                        struct sc_f3_dat  *data)
{
  unsigned int  s, start_1, start_2, length_1, length_2, **a2s;
  int           e_sc, ***sc_up;

  sc_up = data->up_comparative;

  e_sc  = 0;
  a2s   = data->a2s;

  for (s = 0; s < data->n_seq; s++) {
    if (sc_up[s]) {
      start_1   = a2s[s][k] + 1;
      length_1  = a2s[s][l - 1] - a2s[s][k];
      start_2   = a2s[s][i];
      length_2  = a2s[s][i + 1] - a2s[s][i];

      if (length_1 != 0)
        e_sc += sc_up[s][start_1][length_1];

      if (length_2 != 0)
        e_sc += sc_up[s][start_2][length_2];
    }
  }

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_user_split_in_ext_stem1(int              j,
                                 int              k,
                                 int              l,
                                 struct sc_f5_dat *data)
{
  return data->user_cb(1, j, k, l, VRNA_DECOMP_EXT_EXT_STEM1, data->user_data);
}


PRIVATE INLINE int
sc_f3_cb_user_split_in_stem_ext1(int              i,
                                 int              k,
                                 int              l,
                                 struct sc_f3_dat *data)
{
  return data->user_cb(i, data->n, k, l, VRNA_DECOMP_EXT_STEM_EXT1, data->user_data);
}


PRIVATE INLINE int
sc_f5_cb_user_split_in_ext_stem1_comparative(int              j,
                                             int              k,
                                             int              l,
                                             struct sc_f5_dat *data)
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
sc_f3_cb_user_split_in_stem_ext1_comparative(int              i,
                                             int              k,
                                             int              l,
                                             struct sc_f3_dat *data)
{
  unsigned int  s;
  int           e_sc;

  e_sc = 0;

  for (s = 0; s < data->n_seq; s++)
    if (data->user_cb_comparative[s])
      e_sc +=
        data->user_cb_comparative[s](i, data->n, k, l, VRNA_DECOMP_EXT_STEM_EXT1,
                                     data->user_data_comparative[s]);

  return e_sc;
}


PRIVATE INLINE int
sc_f5_cb_user_def_split_in_ext_stem1(int              j,
                                     int              k,
                                     int              l,
                                     struct sc_f5_dat *data)
{
  return sc_f5_cb_split_in_ext_stem1(j, k, l, data) +
         sc_f5_cb_user_split_in_ext_stem1(j, k, l, data);
}


PRIVATE INLINE int
sc_f3_cb_user_def_split_in_stem_ext1(int              i,
                                     int              k,
                                     int              l,
                                     struct sc_f3_dat *data)
{
  return sc_f3_cb_split_in_stem_ext1(i, k, l, data) +
         sc_f3_cb_user_split_in_stem_ext1(i, k, l, data);
}


PRIVATE INLINE int
sc_f5_cb_user_def_split_in_ext_stem1_comparative(int              j,
                                                 int              k,
                                                 int              l,
                                                 struct sc_f5_dat *data)
{
  return sc_f5_cb_split_in_ext_stem1_comparative(j, k, l, data) +
         sc_f5_cb_user_split_in_ext_stem1_comparative(j, k, l, data);
}


PRIVATE INLINE int
sc_f3_cb_user_def_split_in_stem_ext1_comparative(int              i,
                                                 int              k,
                                                 int              l,
                                                 struct sc_f3_dat *data)
{
  return sc_f3_cb_split_in_stem_ext1_comparative(i, k, l, data) +
         sc_f3_cb_user_split_in_stem_ext1_comparative(i, k, l, data);
}
