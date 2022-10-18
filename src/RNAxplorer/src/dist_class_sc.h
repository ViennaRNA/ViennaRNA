#ifndef   _RNAXPLORER_DIST_CLASS_SC_H_
#define   _RNAXPLORER_DIST_CLASS_SC_H_


typedef struct sc_dist_class_s sc_dist_class_t;

typedef double (dist_class_func)(int,
                              int,
                              int,
                              int,
                              unsigned char,
                              int *,
                              sc_dist_class_t *);

typedef void (dist_class_free)(void *);

#include <ViennaRNA/fold_compound.h>

struct sc_dist_class_s {
  unsigned int    ref_num;
  char            **references;
  unsigned int    **ref_bps;
  short           **ref_pts;

  dist_class_func *f;
  void            *f_data;
  dist_class_free *f_free;

  double          kT;
  int             *idx;
};

FLT_OR_DBL
sc_exp_f_dist_class(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char decomp,
                    void          *data);


int
sc_f_dist_class(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char decomp,
                void          *data);


sc_dist_class_t *
sc_dist_class_init(vrna_fold_compound_t *fc);


void
sc_dist_class_destroy(void *data);


void
sc_dist_class_add_ref(sc_dist_class_t *d,
                      const char      *ref_struct);


#endif
