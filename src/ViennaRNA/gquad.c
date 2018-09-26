/*
 * gquad.c
 *
 * Ronny Lorenz 2012
 *
 * ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/gquad.h"

#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif

/**
 *  Use this macro to loop over each G-quadruplex
 *  delimited by a and b within the subsequence [c,d]
 */
#define FOR_EACH_GQUAD(a, b, c, d)  \
  for ((a) = (d) - VRNA_GQUAD_MIN_BOX_SIZE + 1; (a) >= (c); (a)--) \
    for ((b) = (a) + VRNA_GQUAD_MIN_BOX_SIZE - 1; \
         (b) <= MIN2((d), (a) + VRNA_GQUAD_MAX_BOX_SIZE - 1); \
         (b)++)

/**
 *  This macro does almost the same as FOR_EACH_GQUAD() but keeps
 *  the 5' delimiter fixed. 'b' is the 3' delimiter of the gquad,
 *  for gquads within subsequence [a,c] that have 5' delimiter 'a'
 */
#define FOR_EACH_GQUAD_AT(a, b, c)  \
  for ((b) = (a) + VRNA_GQUAD_MIN_BOX_SIZE - 1; \
       (b) <= MIN2((c), (a) + VRNA_GQUAD_MAX_BOX_SIZE - 1); \
       (b)++)


struct gquad_ali_helper {
  short             **S;
  unsigned int      **a2s;
  int               n_seq;
  vrna_param_t      *P;
  vrna_exp_param_t  *pf;
  int               L;
  int               *l;
};

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE INLINE
int *
get_g_islands(short *S);


PRIVATE INLINE
int *
get_g_islands_sub(short *S,
                  int   i,
                  int   j);


/**
 *  IMPORTANT:
 *  If you don't know how to use this function, DONT'T USE IT!
 *
 *  The function pointer this function takes as argument is
 *  used for individual calculations with each g-quadruplex
 *  delimited by [i,j].
 *  The function it points to always receives as first 3 arguments
 *  position i, the stack size L and an array l[3] containing the
 *  individual linker sizes.
 *  The remaining 4 (void *) pointers of the callback function receive
 *  the parameters 'data', 'P', 'aux1' and 'aux2' and thus may be
 *  used to pass whatever data you like to.
 *  As the names of those parameters suggest the convention is that
 *  'data' should be used as a pointer where data is stored into,
 *  e.g the MFE or PF and the 'P' parameter should actually be a
 *  'vrna_param_t *' or 'vrna_exp_param_t *' type.
 *  However, what you actually pass obviously depends on the
 *  function the pointer is pointing to.
 *
 *  Although all of this may look like an overkill, it is found
 *  to be almost as fast as implementing g-quadruplex enumeration
 *  in each individual scenario, i.e. code duplication.
 *  Using this function, however, ensures that all g-quadruplex
 *  enumerations are absolutely identical.
 */
PRIVATE
void
process_gquad_enumeration(int *gg,
                          int i,
                          int j,
                          void (*f)(int, int, int *,
                                    void *, void *, void *, void *),
                          void *data,
                          void *P,
                          void *aux1,
                          void *aux2);


/**
 *  MFE callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe(int   i,
          int   L,
          int   *l,
          void  *data,
          void  *P,
          void  *NA,
          void  *NA2);


PRIVATE
void
gquad_mfe_pos(int   i,
              int   L,
              int   *l,
              void  *data,
              void  *P,
              void  *Lmfe,
              void  *lmfe);


PRIVATE void gquad_mfe_ali_pos(int  i,
                               int  L,
                               int  *l,
                               void *data,
                               void *helper,
                               void *Lmfe,
                               void *lmfe);


PRIVATE
void
gquad_pos_exhaustive(int  i,
                     int  L,
                     int  *l,
                     void *data,
                     void *P,
                     void *Lex,
                     void *lex);


/**
 * Partition function callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_pf(int  i,
         int  L,
         int  *l,
         void *data,
         void *P,
         void *NA,
         void *NA2);


PRIVATE void
gquad_pf_ali(int  i,
             int  L,
             int  *l,
             void *data,
             void *helper,
             void *NA,
             void *NA2);


/**
 * Partition function callback for process_gquad_enumeration()
 * in contrast to gquad_pf() it stores the stack size L and
 * the linker lengths l[3] of the g-quadruplex that dominates
 * the interval [i,j]
 * (FLT_OR_DBL *)data must be 0. on entry
 */
PRIVATE
void
gquad_pf_pos(int  i,
             int  L,
             int  *l,
             void *data,
             void *pf,
             void *Lmax,
             void *lmax);


PRIVATE void
gquad_pf_pos_ali(int  i,
                 int  L,
                 int  *l,
                 void *data,
                 void *helper,
                 void *NA1,
                 void *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe_ali(int   i,
              int   L,
              int   *l,
              void  *data,
              void  *helper,
              void  *NA,
              void  *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 * with seperation of free energy and penalty contribution
 */
PRIVATE
void
gquad_mfe_ali_en(int  i,
                 int  L,
                 int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2);


PRIVATE
void
gquad_interact(int  i,
               int  L,
               int  *l,
               void *data,
               void *pf,
               void *index,
               void *NA2);


PRIVATE
void
gquad_interact_ali(int  i,
                   int  L,
                   int  *l,
                   void *data,
                   void *index,
                   void *helper,
                   void *NA);


PRIVATE
void
gquad_count(int   i,
            int   L,
            int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3);


PRIVATE
void
gquad_count_layers(int  i,
                   int  L,
                   int  *l,
                   void *data,
                   void *NA,
                   void *NA2,
                   void *NA3);


/* other useful static functions */

PRIVATE
int
E_gquad_ali_penalty(int           i,
                    int           L,
                    int           l[3],
                    const short   **S,
                    unsigned int  n_seq,
                    vrna_param_t  *P);


PRIVATE
FLT_OR_DBL
exp_E_gquad_ali_penalty(int               i,
                        int               L,
                        int               l[3],
                        const short       **S,
                        unsigned int      n_seq,
                        vrna_exp_param_t  *P);


PRIVATE void
count_gquad_layer_mismatches(int          i,
                             int          L,
                             int          l[3],
                             const short  **S,
                             unsigned int n_seq,
                             unsigned int mm[2]);


PRIVATE int **
create_L_matrix(short         *S,
                int           start,
                int           maxdist,
                int           n,
                int           **g,
                vrna_param_t  *P);


PRIVATE int **
create_aliL_matrix(int          start,
                   int          maxdist,
                   int          n,
                   int          **g,
                   short        *S_cons,
                   short        **S,
                   unsigned int **a2s,
                   int          n_seq,
                   vrna_param_t *P);


PRIVATE void gquad_mfe_ali_pos(int  i,
                               int  L,
                               int  *l,
                               void *data,
                               void *P,
                               void *Lmfe,
                               void *lmfe);


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */

/********************************
 * Here are the G-quadruplex energy
 * contribution functions
 *********************************/
PUBLIC int
E_gquad(int           L,
        int           l[3],
        vrna_param_t  *P)
{
  int i, c = INF;

  for (i = 0; i < 3; i++) {
    if (l[i] > VRNA_GQUAD_MAX_LINKER_LENGTH)
      return c;

    if (l[i] < VRNA_GQUAD_MIN_LINKER_LENGTH)
      return c;
  }
  if (L > VRNA_GQUAD_MAX_STACK_SIZE)
    return c;

  if (L < VRNA_GQUAD_MIN_STACK_SIZE)
    return c;

  gquad_mfe(0, L, l,
            (void *)(&c),
            (void *)P,
            NULL,
            NULL);
  return c;
}


PUBLIC FLT_OR_DBL
exp_E_gquad(int               L,
            int               l[3],
            vrna_exp_param_t  *pf)
{
  int         i;
  FLT_OR_DBL  q = 0.;

  for (i = 0; i < 3; i++) {
    if (l[i] > VRNA_GQUAD_MAX_LINKER_LENGTH)
      return q;

    if (l[i] < VRNA_GQUAD_MIN_LINKER_LENGTH)
      return q;
  }
  if (L > VRNA_GQUAD_MAX_STACK_SIZE)
    return q;

  if (L < VRNA_GQUAD_MIN_STACK_SIZE)
    return q;

  gquad_pf(0, L, l,
           (void *)(&q),
           (void *)pf,
           NULL,
           NULL);
  return q;
}


PUBLIC FLT_OR_DBL
exp_E_gquad_ali(int               i,
                int               L,
                int               l[3],
                short             **S,
                unsigned int      **a2s,
                int               n_seq,
                vrna_exp_param_t  *pf)
{
  int         s;
  FLT_OR_DBL  q = 0.;

  for (s = 0; s < 3; s++) {
    if (l[s] > VRNA_GQUAD_MAX_LINKER_LENGTH)
      return q;

    if (l[s] < VRNA_GQUAD_MIN_LINKER_LENGTH)
      return q;
  }
  if (L > VRNA_GQUAD_MAX_STACK_SIZE)
    return q;

  if (L < VRNA_GQUAD_MIN_STACK_SIZE)
    return q;

  struct gquad_ali_helper gq_help;

  gq_help.S     = S;
  gq_help.a2s   = a2s;
  gq_help.n_seq = n_seq;
  gq_help.pf    = pf;

  gquad_pf_ali(i, L, l,
               (void *)(&q),
               (void *)&gq_help,
               NULL,
               NULL);
  return q;
}


PUBLIC void
E_gquad_ali_en(int          i,
               int          L,
               int          l[3],
               const short  **S,
               unsigned int **a2s,
               unsigned int n_seq,
               vrna_param_t *P,
               int          en[2])
{
  unsigned int  s;
  int           ee, ee2, u1, u2, u3;

  en[0] = en[1] = INF;

  /* check if the quadruplex obeys the canonical form */
  for (s = 0; s < 3; s++) {
    if (l[s] > VRNA_GQUAD_MAX_LINKER_LENGTH)
      return;

    if (l[s] < VRNA_GQUAD_MIN_LINKER_LENGTH)
      return;
  }
  if (L > VRNA_GQUAD_MAX_STACK_SIZE)
    return;

  if (L < VRNA_GQUAD_MIN_STACK_SIZE)
    return;

  /* compute actual quadruplex contribution for subalignment */
  ee = 0;

  for (s = 0; s < n_seq; s++) {
    u1  = a2s[s][i + L + l[0] - 1] - a2s[s][i + L - 1];
    u2  = a2s[s][i + 2 * L + l[0] + l[1] - 1] - a2s[s][i + 2 * L + l[0] - 1];
    u3  = a2s[s][i + 3 * L + l[0] + l[1] + l[2] - 1] - a2s[s][i + 3 * L + l[0] + l[1] - 1];
    ee  += P->gquad[L][u1 + u2 + u3];
  }

  /* get penalty from incompatible sequences in alignment */
  ee2 = E_gquad_ali_penalty(i, L, l, S, n_seq, P);

  /* assign return values */
  if (ee2 != INF) {
    en[0] = ee;
    en[1] = ee2;
  }
}


/********************************
 * Now, the triangular matrix
 * generators for the G-quadruplex
 * contributions are following
 *********************************/
PUBLIC int *
get_gquad_matrix(short        *S,
                 vrna_param_t *P)
{
  int n, size, i, j, *gg, *my_index, *data;

  n         = S[0];
  my_index  = vrna_idx_col_wise(n);
  gg        = get_g_islands(S);
  size      = (n * (n + 1)) / 2 + 2;
  data      = (int *)vrna_alloc(sizeof(int) * size);

  /* prefill the upper triangular matrix with INF */
  for (i = 0; i < size; i++)
    data[i] = INF;

  FOR_EACH_GQUAD(i, j, 1, n){
    process_gquad_enumeration(gg, i, j,
                              &gquad_mfe,
                              (void *)(&(data[my_index[j] + i])),
                              (void *)P,
                              NULL,
                              NULL);
  }

  free(my_index);
  free(gg);
  return data;
}


PUBLIC FLT_OR_DBL *
get_gquad_pf_matrix(short             *S,
                    FLT_OR_DBL        *scale,
                    vrna_exp_param_t  *pf)
{
  int         n, size, *gg, i, j, *my_index;
  FLT_OR_DBL  *data;


  n         = S[0];
  size      = (n * (n + 1)) / 2 + 2;
  data      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  gg        = get_g_islands(S);
  my_index  = vrna_idx_row_wise(n);

  FOR_EACH_GQUAD(i, j, 1, n){
    process_gquad_enumeration(gg, i, j,
                              &gquad_pf,
                              (void *)(&(data[my_index[i] - j])),
                              (void *)pf,
                              NULL,
                              NULL);
    data[my_index[i] - j] *= scale[j - i + 1];
  }

  free(my_index);
  free(gg);
  return data;
}


PUBLIC FLT_OR_DBL *
get_gquad_pf_matrix_comparative(short             *S_cons,
                                short             **S,
                                unsigned int      **a2s,
                                FLT_OR_DBL        *scale,
                                unsigned int      n_seq,
                                vrna_exp_param_t  *pf)
{
  int                     n, size, *gg, i, j, *my_index;
  FLT_OR_DBL              *data;
  struct gquad_ali_helper gq_help;


  n             = S[0][0];
  size          = (n * (n + 1)) / 2 + 2;
  data          = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  gg            = get_g_islands(S_cons);
  my_index      = vrna_idx_row_wise(n);
  gq_help.S     = S;
  gq_help.a2s   = a2s;
  gq_help.n_seq = n_seq;
  gq_help.pf    = pf;

  FOR_EACH_GQUAD(i, j, 1, n){
    process_gquad_enumeration(gg, i, j,
                              &gquad_pf_ali,
                              (void *)(&(data[my_index[i] - j])),
                              (void *)&gq_help,
                              NULL,
                              NULL);
    data[my_index[i] - j] *= scale[j - i + 1];
  }

  free(my_index);
  free(gg);
  return data;
}


PUBLIC int *
get_gquad_ali_matrix(short        *S_cons,
                     short        **S,
                     unsigned int **a2s,
                     int          n_seq,
                     vrna_param_t *P)
{
  int                     n, size, *data, *gg;
  int                     i, j, *my_index;
  struct gquad_ali_helper gq_help;

  n         = S[0][0];
  size      = (n * (n + 1)) / 2 + 2;
  data      = (int *)vrna_alloc(sizeof(int) * size);
  gg        = get_g_islands(S_cons);
  my_index  = vrna_idx_col_wise(n);

  gq_help.S     = S;
  gq_help.a2s   = a2s;
  gq_help.n_seq = n_seq;
  gq_help.P     = P;

  /* prefill the upper triangular matrix with INF */
  for (i = 0; i < size; i++)
    data[i] = INF;

  FOR_EACH_GQUAD(i, j, 1, n){
    process_gquad_enumeration(gg, i, j,
                              &gquad_mfe_ali,
                              (void *)(&(data[my_index[j] + i])),
                              (void *)&gq_help,
                              NULL,
                              NULL);
  }

  free(my_index);
  free(gg);
  return data;
}


PUBLIC int **
get_gquad_L_matrix(short        *S,
                   int          start,
                   int          maxdist,
                   int          n,
                   int          **g,
                   vrna_param_t *P)
{
  return create_L_matrix(S, start, maxdist, n, g, P);
}


PUBLIC void
vrna_gquad_mx_local_update(vrna_fold_compound_t *vc,
                           int                  start)
{
  if (vc->type == VRNA_FC_TYPE_COMPARATIVE) {
    vc->matrices->ggg_local = create_aliL_matrix(
      start,
      vc->window_size,
      vc->length,
      vc->matrices->ggg_local,
      vc->S_cons,
      vc->S,
      vc->a2s,
      vc->n_seq,
      vc->params);
  } else {
    vc->matrices->ggg_local = create_L_matrix(
      vc->sequence_encoding,
      start,
      vc->window_size,
      vc->length,
      vc->matrices->ggg_local,
      vc->params);
  }
}


PRIVATE int **
create_L_matrix(short         *S,
                int           start,
                int           maxdist,
                int           n,
                int           **g,
                vrna_param_t  *P)
{
  int **data;
  int i, j, k, *gg, p, q;

  p   = MAX2(1, start);
  q   = MIN2(n, start + maxdist + 4);
  gg  = get_g_islands_sub(S, p, q);

  if (g) {
    /* we just update the gquadruplex contribution for the current
     * start and rotate the rest */
    data = g;
    /* we re-use the memory allocated previously */
    data[start]               = data[start + maxdist + 5];
    data[start + maxdist + 5] = NULL;

    /* prefill with INF */
    for (i = 0; i < maxdist + 5; i++)
      data[start][i] = INF;

    /*  now we compute contributions for all gquads with 5' delimiter at
     *  position 'start'
     */
    FOR_EACH_GQUAD_AT(start, j, start + maxdist + 4){
      process_gquad_enumeration(gg, start, j,
                                &gquad_mfe,
                                (void *)(&(data[start][j - start])),
                                (void *)P,
                                NULL,
                                NULL);
    }
  } else {
    /* create a new matrix from scratch since this is the first
     * call to this function */

    /* allocate memory and prefill with INF */
    data = (int **)vrna_alloc(sizeof(int *) * (n + 1));
    for (k = n; (k > n - maxdist - 5) && (k >= 0); k--) {
      data[k] = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
      for (i = 0; i < maxdist + 5; i++)
        data[k][i] = INF;
    }

    /* compute all contributions for the gquads in this interval */
    FOR_EACH_GQUAD(i, j, MAX2(1, n - maxdist - 4), n){
      process_gquad_enumeration(gg, i, j,
                                &gquad_mfe,
                                (void *)(&(data[i][j - i])),
                                (void *)P,
                                NULL,
                                NULL);
    }
  }

  gg += p - 1;
  free(gg);
  return data;
}


PRIVATE int **
create_aliL_matrix(int          start,
                   int          maxdist,
                   int          n,
                   int          **g,
                   short        *S_cons,
                   short        **S,
                   unsigned int **a2s,
                   int          n_seq,
                   vrna_param_t *P)
{
  int **data;
  int i, j, k, *gg, p, q;

  p   = MAX2(1, start);
  q   = MIN2(n, start + maxdist + 4);
  gg  = get_g_islands_sub(S_cons, p, q);

  struct gquad_ali_helper gq_help;

  gq_help.S     = S;
  gq_help.a2s   = a2s;
  gq_help.n_seq = n_seq;
  gq_help.P     = P;

  if (g) {
    /* we just update the gquadruplex contribution for the current
     * start and rotate the rest */
    data = g;
    /* we re-use the memory allocated previously */
    data[start]               = data[start + maxdist + 5];
    data[start + maxdist + 5] = NULL;

    /* prefill with INF */
    for (i = 0; i < maxdist + 5; i++)
      data[start][i] = INF;

    /*  now we compute contributions for all gquads with 5' delimiter at
     *  position 'start'
     */
    FOR_EACH_GQUAD_AT(start, j, start + maxdist + 4){
      process_gquad_enumeration(gg, start, j,
                                &gquad_mfe_ali,
                                (void *)(&(data[start][j - start])),
                                (void *)&gq_help,
                                NULL,
                                NULL);
    }
  } else {
    /* create a new matrix from scratch since this is the first
     * call to this function */

    /* allocate memory and prefill with INF */
    data = (int **)vrna_alloc(sizeof(int *) * (n + 1));
    for (k = n; (k > n - maxdist - 5) && (k >= 0); k--) {
      data[k] = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
      for (i = 0; i < maxdist + 5; i++)
        data[k][i] = INF;
    }

    /* compute all contributions for the gquads in this interval */
    FOR_EACH_GQUAD(i, j, MAX2(1, n - maxdist - 4), n){
      process_gquad_enumeration(gg, i, j,
                                &gquad_mfe_ali,
                                (void *)(&(data[i][j - i])),
                                (void *)&gq_help,
                                NULL,
                                NULL);
    }
  }

  gg += p - 1;
  free(gg);
  return data;
}


PUBLIC plist *
get_plist_gquad_from_db(const char  *structure,
                        float       pr)
{
  int   x, size, actual_size, L, n, ge, ee, gb, l[3];
  plist *pl;

  actual_size = 0;
  ge          = 0;
  n           = 2;
  size        = strlen(structure);
  pl          = (plist *)vrna_alloc(n * size * sizeof(plist));

  while ((ee = parse_gquad(structure + ge, &L, l)) > 0) {
    ge  += ee;
    gb  = ge - L * 4 - l[0] - l[1] - l[2] + 1;
    /* add pseudo-base pair encloding gquad */
    for (x = 0; x < L; x++) {
      if (actual_size >= n * size - 5) {
        n   *= 2;
        pl  = (plist *)vrna_realloc(pl, n * size * sizeof(plist));
      }

      pl[actual_size].i       = gb + x;
      pl[actual_size].j       = ge + x - L + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = 0;

      pl[actual_size].i       = gb + x;
      pl[actual_size].j       = gb + x + l[0] + L;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = 0;

      pl[actual_size].i       = gb + x + l[0] + L;
      pl[actual_size].j       = ge + x - 2 * L - l[2] + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = 0;

      pl[actual_size].i       = ge + x - 2 * L - l[2] + 1;
      pl[actual_size].j       = ge + x - L + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = 0;
    }
  }

  pl[actual_size].i   = pl[actual_size].j = 0;
  pl[actual_size++].p = 0;
  pl                  = (plist *)vrna_realloc(pl, actual_size * sizeof(plist));
  return pl;
}


PUBLIC void
get_gquad_pattern_mfe(short         *S,
                      int           i,
                      int           j,
                      vrna_param_t  *P,
                      int           *L,
                      int           l[3])
{
  int *gg = get_g_islands_sub(S, i, j);
  int c   = INF;

  process_gquad_enumeration(gg, i, j,
                            &gquad_mfe_pos,
                            (void *)(&c),
                            (void *)P,
                            (void *)L,
                            (void *)l);

  gg += i - 1;
  free(gg);
}


PUBLIC void
get_gquad_pattern_mfe_ali(short         **S,
                          unsigned int  **a2s,
                          short         *S_cons,
                          int           n_seq,
                          int           i,
                          int           j,
                          vrna_param_t  *P,
                          int           *L,
                          int           l[3])
{
  int                     mfe, *gg;
  struct gquad_ali_helper gq_help;

  gg  = get_g_islands_sub(S_cons, i, j);
  mfe = INF;

  gq_help.S     = S;
  gq_help.a2s   = a2s;
  gq_help.n_seq = n_seq;
  gq_help.P     = P;

  process_gquad_enumeration(gg, i, j,
                            &gquad_mfe_ali_pos,
                            (void *)(&mfe),
                            (void *)&gq_help,
                            (void *)L,
                            (void *)l);

  gg += i - 1;
  free(gg);
}


PUBLIC void
get_gquad_pattern_exhaustive(short        *S,
                             int          i,
                             int          j,
                             vrna_param_t *P,
                             int          *L,
                             int          *l,
                             int          threshold)
{
  int *gg = get_g_islands_sub(S, i, j);

  process_gquad_enumeration(gg, i, j,
                            &gquad_pos_exhaustive,
                            (void *)(&threshold),
                            (void *)P,
                            (void *)L,
                            (void *)l);

  gg += i - 1;
  free(gg);
}


PUBLIC void
get_gquad_pattern_pf(short            *S,
                     int              i,
                     int              j,
                     vrna_exp_param_t *pf,
                     int              *L,
                     int              l[3])
{
  int         *gg = get_g_islands_sub(S, i, j);
  FLT_OR_DBL  q   = 0.;

  process_gquad_enumeration(gg, i, j,
                            &gquad_pf_pos,
                            (void *)(&q),
                            (void *)pf,
                            (void *)L,
                            (void *)l);

  gg += i - 1;
  free(gg);
}


PUBLIC void
vrna_get_gquad_pattern_pf(vrna_fold_compound_t  *fc,
                          int                   i,
                          int                   j,
                          int                   *L,
                          int                   l[3])
{
  short             *S  = fc->type == VRNA_FC_TYPE_SINGLE ? fc->sequence_encoding2 : fc->S_cons;
  int               *gg = get_g_islands_sub(S, i, j);
  FLT_OR_DBL        q   = 0.;
  vrna_exp_param_t  *pf = fc->exp_params;

  if (fc->type == VRNA_FC_TYPE_SINGLE) {
    process_gquad_enumeration(gg, i, j,
                              &gquad_pf_pos,
                              (void *)(&q),
                              (void *)pf,
                              (void *)L,
                              (void *)l);
  } else {
    struct gquad_ali_helper gq_help;
    gq_help.S     = fc->S;
    gq_help.a2s   = fc->a2s;
    gq_help.n_seq = fc->n_seq;
    gq_help.pf    = pf;
    gq_help.L     = (int)(*L);
    gq_help.l     = (int *)l;
    process_gquad_enumeration(gg, i, j,
                              &gquad_pf_pos_ali,
                              (void *)(&q),
                              (void *)&gq_help,
                              NULL,
                              NULL);
    *L = gq_help.L;
  }

  gg += i - 1;
  free(gg);
}


PUBLIC plist *
get_plist_gquad_from_pr(short             *S,
                        int               gi,
                        int               gj,
                        FLT_OR_DBL        *G,
                        FLT_OR_DBL        *probs,
                        FLT_OR_DBL        *scale,
                        vrna_exp_param_t  *pf)
{
  int L, l[3];

  return get_plist_gquad_from_pr_max(S, gi, gj, G, probs, scale, &L, l, pf);
}


PUBLIC plist *
vrna_get_plist_gquad_from_pr(vrna_fold_compound_t *fc,
                             int                  gi,
                             int                  gj)
{
  int L, l[3];

  return vrna_get_plist_gquad_from_pr_max(fc, gi, gj, &L, l);
}


PUBLIC plist *
get_plist_gquad_from_pr_max(short             *S,
                            int               gi,
                            int               gj,
                            FLT_OR_DBL        *G,
                            FLT_OR_DBL        *probs,
                            FLT_OR_DBL        *scale,
                            int               *Lmax,
                            int               lmax[3],
                            vrna_exp_param_t  *pf)
{
  int         n, size, *gg, counter, i, j, *my_index;
  FLT_OR_DBL  pp, *tempprobs;
  plist       *pl;

  n         = S[0];
  size      = (n * (n + 1)) / 2 + 2;
  tempprobs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  pl        = (plist *)vrna_alloc((S[0] * S[0]) * sizeof(plist));
  gg        = get_g_islands_sub(S, gi, gj);
  counter   = 0;
  my_index  = vrna_idx_row_wise(n);

  process_gquad_enumeration(gg, gi, gj,
                            &gquad_interact,
                            (void *)tempprobs,
                            (void *)pf,
                            (void *)my_index,
                            NULL);

  pp = 0.;
  process_gquad_enumeration(gg, gi, gj,
                            &gquad_pf_pos,
                            (void *)(&pp),
                            (void *)pf,
                            (void *)Lmax,
                            (void *)lmax);

  pp = probs[my_index[gi] - gj] * scale[gj - gi + 1] / G[my_index[gi] - gj];
  for (i = gi; i < gj; i++) {
    for (j = i; j <= gj; j++) {
      if (tempprobs[my_index[i] - j] > 0.) {
        pl[counter].i   = i;
        pl[counter].j   = j;
        pl[counter++].p = pp * tempprobs[my_index[i] - j];
      }
    }
  }
  pl[counter].i   = pl[counter].j = 0;
  pl[counter++].p = 0.;
  /* shrink memory to actual size needed */
  pl = (plist *)vrna_realloc(pl, counter * sizeof(plist));

  gg += gi - 1;
  free(gg);
  free(my_index);
  free(tempprobs);
  return pl;
}


PUBLIC plist *
vrna_get_plist_gquad_from_pr_max(vrna_fold_compound_t *fc,
                                 int                  gi,
                                 int                  gj,
                                 int                  *Lmax,
                                 int                  lmax[3])
{
  short             *S;
  int               n, size, *gg, counter, i, j, *my_index;
  FLT_OR_DBL        pp, *tempprobs, *G, *probs, *scale;
  plist             *pl;
  vrna_exp_param_t  *pf;

  n         = (int)fc->length;
  pf        = fc->exp_params;
  G         = fc->exp_matrices->G;
  probs     = fc->exp_matrices->probs;
  scale     = fc->exp_matrices->scale;
  S         = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : fc->S_cons;
  size      = (n * (n + 1)) / 2 + 2;
  tempprobs = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  pl        = (plist *)vrna_alloc((n * n) * sizeof(plist));
  gg        = get_g_islands_sub(S, gi, gj);
  counter   = 0;
  my_index  = vrna_idx_row_wise(n);

  pp = 0.;

  if (fc->type == VRNA_FC_TYPE_SINGLE) {
    process_gquad_enumeration(gg, gi, gj,
                              &gquad_interact,
                              (void *)tempprobs,
                              (void *)pf,
                              (void *)my_index,
                              NULL);

    process_gquad_enumeration(gg, gi, gj,
                              &gquad_pf_pos,
                              (void *)(&pp),
                              (void *)pf,
                              (void *)Lmax,
                              (void *)lmax);
  } else {
    struct gquad_ali_helper gq_help;
    gq_help.S     = fc->S;
    gq_help.a2s   = fc->a2s;
    gq_help.n_seq = fc->n_seq;
    gq_help.pf    = pf;
    gq_help.L     = *Lmax;
    gq_help.l     = &(lmax[0]);

    process_gquad_enumeration(gg, gi, gj,
                              &gquad_interact_ali,
                              (void *)tempprobs,
                              (void *)my_index,
                              (void *)&(gq_help),
                              NULL);

    process_gquad_enumeration(gg, gi, gj,
                              &gquad_pf_pos_ali,
                              (void *)(&pp),
                              (void *)&gq_help,
                              NULL,
                              NULL);
    *Lmax = gq_help.L;
  }

  pp = probs[my_index[gi] - gj] * scale[gj - gi + 1] / G[my_index[gi] - gj];
  for (i = gi; i < gj; i++) {
    for (j = i; j <= gj; j++) {
      if (tempprobs[my_index[i] - j] > 0.) {
        pl[counter].i   = i;
        pl[counter].j   = j;
        pl[counter++].p = pp * tempprobs[my_index[i] - j];
      }
    }
  }
  pl[counter].i   = pl[counter].j = 0;
  pl[counter++].p = 0.;
  /* shrink memory to actual size needed */
  pl = (plist *)vrna_realloc(pl, counter * sizeof(plist));

  gg += gi - 1;
  free(gg);
  free(my_index);
  free(tempprobs);
  return pl;
}


PUBLIC int
get_gquad_count(short *S,
                int   i,
                int   j)
{
  int *gg = get_g_islands_sub(S, i, j);
  int p, q, counter = 0;

  FOR_EACH_GQUAD(p, q, i, j)
  process_gquad_enumeration(gg, p, q,
                            &gquad_count,
                            (void *)(&counter),
                            NULL,
                            NULL,
                            NULL);

  gg += i - 1;
  free(gg);
  return counter;
}


PUBLIC int
get_gquad_layer_count(short *S,
                      int   i,
                      int   j)
{
  int *gg = get_g_islands_sub(S, i, j);
  int p, q, counter = 0;

  FOR_EACH_GQUAD(p, q, i, j)
  process_gquad_enumeration(gg, p, q,
                            &gquad_count_layers,
                            (void *)(&counter),
                            NULL,
                            NULL,
                            NULL);

  gg += i - 1;
  free(gg);
  return counter;
}


PUBLIC int
parse_gquad(const char  *struc,
            int         *L,
            int         l[3])
{
  int i, il, start, end, len;

  for (i = 0; struc[i] && struc[i] != '+'; i++);
  if (struc[i] == '+') {
    /* start of gquad */
    for (il = 0; il <= 3; il++) {
      start = i; /* pos of first '+' */
      while (struc[++i] == '+')
        if ((il) && (i - start == *L))
          break;

      end = i;
      len = end - start;
      if (il == 0)
        *L = len;
      else if (len != *L)
        vrna_message_error("unequal stack lengths in gquad");

      if (il == 3)
        break;

      while (struc[++i] == '.'); /* linker */
      l[il] = i - end;
      if (struc[i] != '+')
        vrna_message_error("illegal character in gquad linker region");
    }
  } else {
    return 0;
  }

  /* printf("gquad at %d %d %d %d %d\n", end, *L, l[0], l[1], l[2]); */
  return end;
}


/*
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */
PRIVATE int
E_gquad_ali_penalty(int           i,
                    int           L,
                    int           l[3],
                    const short   **S,
                    unsigned int  n_seq,
                    vrna_param_t  *P)
{
  unsigned int mm[2];

  count_gquad_layer_mismatches(i, L, l, S, n_seq, mm);

  if (mm[1] > P->gquadLayerMismatchMax)
    return INF;
  else
    return P->gquadLayerMismatch * mm[0];
}


PRIVATE FLT_OR_DBL
exp_E_gquad_ali_penalty(int               i,
                        int               L,
                        int               l[3],
                        const short       **S,
                        unsigned int      n_seq,
                        vrna_exp_param_t  *pf)
{
  unsigned int mm[2];

  count_gquad_layer_mismatches(i, L, l, S, n_seq, mm);

  if (mm[1] > pf->gquadLayerMismatchMax)
    return (FLT_OR_DBL)0.;
  else
    return (FLT_OR_DBL)pow(pf->expgquadLayerMismatch, (double)mm[0]);
}


PRIVATE void
count_gquad_layer_mismatches(int          i,
                             int          L,
                             int          l[3],
                             const short  **S,
                             unsigned int n_seq,
                             unsigned int mm[2])
{
  unsigned int  s;
  int           cnt;

  mm[0] = mm[1] = 0;

  /* check for compatibility in the alignment */
  for (s = 0; s < n_seq; s++) {
    unsigned int  ld        = 0; /* !=0 if layer destruction was detected */
    unsigned int  mismatch  = 0;

    /* check bottom layer */
    if (S[s][i] != 3)
      ld |= 1U;

    if (S[s][i + L + l[0]] != 3)
      ld |= 2U;

    if (S[s][i + 2 * L + l[0] + l[1]] != 3)
      ld |= 4U;

    if (S[s][i + 3 * L + l[0] + l[1] + l[2]] != 3)
      ld |= 8U;

    /* add 1x penalty for missing bottom layer */
    if (ld)
      mismatch++;

    /* check top layer */
    ld = 0;
    if (S[s][i + L - 1] != 3)
      ld |= 1U;

    if (S[s][i + 2 * L + l[0] - 1] != 3)
      ld |= 2U;

    if (S[s][i + 3 * L + l[0] + l[1] - 1] != 3)
      ld |= 4U;

    if (S[s][i + 4 * L + l[0] + l[1] + l[2] - 1] != 3)
      ld |= 8U;

    /* add 1x penalty for missing top layer */
    if (ld)
      mismatch++;

    /* check inner layers */
    ld = 0;
    for (cnt = 1; cnt < L - 1; cnt++) {
      if (S[s][i + cnt] != 3)
        ld |= 1U;

      if (S[s][i + L + l[0] + cnt] != 3)
        ld |= 2U;

      if (S[s][i + 2 * L + l[0] + l[1] + cnt] != 3)
        ld |= 4U;

      if (S[s][i + 3 * L + l[0] + l[1] + l[2] + cnt] != 3)
        ld |= 8U;

      /* add 2x penalty for missing inner layer */
      if (ld)
        mismatch += 2;
    }

    mm[0] += mismatch;

    if (mismatch >= 2 * (L - 1))
      mm[1]++;
  }
}


PRIVATE void
gquad_mfe(int   i,
          int   L,
          int   *l,
          void  *data,
          void  *P,
          void  *NA,
          void  *NA2)
{
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc < *((int *)data))
    *((int *)data) = cc;
}


PRIVATE void
gquad_mfe_pos(int   i,
              int   L,
              int   *l,
              void  *data,
              void  *P,
              void  *Lmfe,
              void  *lmfe)
{
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc < *((int *)data)) {
    *((int *)data)        = cc;
    *((int *)Lmfe)        = L;
    *((int *)lmfe)        = l[0];
    *(((int *)lmfe) + 1)  = l[1];
    *(((int *)lmfe) + 2)  = l[2];
  }
}


PRIVATE void
gquad_mfe_ali_pos(int   i,
                  int   L,
                  int   *l,
                  void  *data,
                  void  *helper,
                  void  *Lmfe,
                  void  *lmfe)
{
  int cc = INF;

  gquad_mfe_ali(i, L, l, (void *)&cc, helper, NULL, NULL);

  if (cc < *((int *)data)) {
    *((int *)data)        = cc;
    *((int *)Lmfe)        = L;
    *((int *)lmfe)        = l[0];
    *(((int *)lmfe) + 1)  = l[1];
    *(((int *)lmfe) + 2)  = l[2];
  }
}


PRIVATE
void
gquad_pos_exhaustive(int  i,
                     int  L,
                     int  *l,
                     void *data,
                     void *P,
                     void *Lex,
                     void *lex)
{
  int cnt;
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc <= *((int *)data)) {
    /*  since Lex is an array of L values and lex an
     *  array of l triples we need to find out where
     *  the current gquad position is to be stored...
     * the below implementation might be slow but we
     * still use it for now
     */
    for (cnt = 0; ((int *)Lex)[cnt] != -1; cnt++);

    *((int *)Lex + cnt)             = L;
    *((int *)Lex + cnt + 1)         = -1;
    *(((int *)lex) + (3 * cnt) + 0) = l[0];
    *(((int *)lex) + (3 * cnt) + 1) = l[1];
    *(((int *)lex) + (3 * cnt) + 2) = l[2];
  }
}


PRIVATE
void
gquad_count(int   i,
            int   L,
            int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3)
{
  *((int *)data) += 1;
}


PRIVATE
void
gquad_count_layers(int  i,
                   int  L,
                   int  *l,
                   void *data,
                   void *NA,
                   void *NA2,
                   void *NA3)
{
  *((int *)data) += L;
}


PRIVATE void
gquad_pf(int  i,
         int  L,
         int  *l,
         void *data,
         void *pf,
         void *NA,
         void *NA2)
{
  *((FLT_OR_DBL *)data) += ((vrna_exp_param_t *)pf)->expgquad[L][l[0] + l[1] + l[2]];
}


PRIVATE void
gquad_pf_ali(int  i,
             int  L,
             int  *l,
             void *data,
             void *helper,
             void *NA,
             void *NA2)
{
  short                   **S;
  unsigned int            **a2s;
  int                     u1, u2, u3, s, n_seq;
  FLT_OR_DBL              penalty;
  vrna_exp_param_t        *pf;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n_seq   = gq_help->n_seq;
  pf      = gq_help->pf;
  penalty = exp_E_gquad_ali_penalty(i, L, l, (const short **)S, (unsigned int)n_seq, pf);

  if (penalty != 0.) {
    double q = 1.;
    for (s = 0; s < n_seq; s++) {
      u1  = a2s[s][i + L + l[0] - 1] - a2s[s][i + L - 1];
      u2  = a2s[s][i + 2 * L + l[0] + l[1] - 1] - a2s[s][i + 2 * L + l[0] - 1];
      u3  = a2s[s][i + 3 * L + l[0] + l[1] + l[2] - 1] - a2s[s][i + 3 * L + l[0] + l[1] - 1];
      q   *= pf->expgquad[L][u1 + u2 + u3];
    }
    *((FLT_OR_DBL *)data) += q * penalty;
  }
}


PRIVATE void
gquad_pf_pos(int  i,
             int  L,
             int  *l,
             void *data,
             void *pf,
             void *Lmax,
             void *lmax)
{
  FLT_OR_DBL gq = 0.;

  gquad_pf(i, L, l, (void *)&gq, pf, NULL, NULL);

  if (gq > *((FLT_OR_DBL *)data)) {
    *((FLT_OR_DBL *)data) = gq;
    *((int *)Lmax)        = L;
    *((int *)lmax)        = l[0];
    *(((int *)lmax) + 1)  = l[1];
    *(((int *)lmax) + 2)  = l[2];
  }
}


PRIVATE void
gquad_pf_pos_ali(int  i,
                 int  L,
                 int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2)
{
  FLT_OR_DBL              gq        = 0.;
  struct gquad_ali_helper *gq_help  = (struct gquad_ali_helper *)helper;

  gquad_pf_ali(i, L, l, (void *)&gq, helper, NULL, NULL);

  if (gq > *((FLT_OR_DBL *)data)) {
    *((FLT_OR_DBL *)data) = gq;
    gq_help->L            = L;
    gq_help->l[0]         = l[0];
    gq_help->l[1]         = l[1];
    gq_help->l[2]         = l[2];
  }
}


PRIVATE void
gquad_mfe_ali(int   i,
              int   L,
              int   *l,
              void  *data,
              void  *helper,
              void  *NA,
              void  *NA2)
{
  int j, en[2], cc;

  en[0] = en[1] = INF;

  for (j = 0; j < 3; j++) {
    if (l[j] > VRNA_GQUAD_MAX_LINKER_LENGTH)
      return;

    if (l[j] < VRNA_GQUAD_MIN_LINKER_LENGTH)
      return;
  }
  if (L > VRNA_GQUAD_MAX_STACK_SIZE)
    return;

  if (L < VRNA_GQUAD_MIN_STACK_SIZE)
    return;

  gquad_mfe_ali_en(i, L, l, (void *)(&(en[0])), helper, NULL, NULL);
  if (en[1] != INF) {
    cc = en[0] + en[1];
    if (cc < *((int *)data))
      *((int *)data) = cc;
  }
}


PRIVATE void
gquad_mfe_ali_en(int  i,
                 int  L,
                 int  *l,
                 void *data,
                 void *helper,
                 void *NA,
                 void *NA2)
{
  short                   **S;
  unsigned int            **a2s;
  int                     en[2], cc, dd, s, n_seq, u1, u2, u3;
  vrna_param_t            *P;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n_seq   = gq_help->n_seq;
  P       = gq_help->P;

  for (en[0] = s = 0; s < n_seq; s++) {
    u1    = a2s[s][i + L + l[0] - 1] - a2s[s][i + L - 1];
    u2    = a2s[s][i + 2 * L + l[0] + l[1] - 1] - a2s[s][i + 2 * L + l[0] - 1];
    u3    = a2s[s][i + 3 * L + l[0] + l[1] + l[2] - 1] - a2s[s][i + 3 * L + l[0] + l[1] - 1];
    en[0] += P->gquad[L][u1 + u2 + u3];
  }
  en[1] = E_gquad_ali_penalty(i, L, l, (const short **)S, (unsigned int)n_seq, P);

  if (en[1] != INF) {
    cc  = en[0] + en[1];
    dd  = ((int *)data)[0] + ((int *)data)[1];
    if (cc < dd) {
      ((int *)data)[0]  = en[0];
      ((int *)data)[1]  = en[1];
    }
  }
}


PRIVATE void
gquad_interact(int  i,
               int  L,
               int  *l,
               void *data,
               void *pf,
               void *index,
               void *NA2)
{
  int         x, *idx;
  FLT_OR_DBL  gq, *pp;

  idx = (int *)index;
  pp  = (FLT_OR_DBL *)data;
  gq  = exp_E_gquad(L, l, (vrna_exp_param_t *)pf);

  for (x = 0; x < L; x++) {
    pp[idx[i + x] - (i + x + 3 * L + l[0] + l[1] + l[2])]                       += gq;
    pp[idx[i + x] - (i + x + L + l[0])]                                         += gq;
    pp[idx[i + x + L + l[0]] - (i + x + 2 * L + l[0] + l[1])]                   += gq;
    pp[idx[i + x + 2 * L + l[0] + l[1]] - (i + x + 3 * L + l[0] + l[1] + l[2])] += gq;
  }
}


PRIVATE void
gquad_interact_ali(int  i,
                   int  L,
                   int  *l,
                   void *data,
                   void *index,
                   void *helper,
                   void *NA)
{
  int         x, *idx, bad;
  FLT_OR_DBL  gq, *pp;

  idx = (int *)index;
  pp  = (FLT_OR_DBL *)data;
  bad = 0;

  for (x = 0; x < 3; x++) {
    if (l[x] > VRNA_GQUAD_MAX_LINKER_LENGTH) {
      bad = 1;
      break;
    }

    if (l[x] < VRNA_GQUAD_MIN_LINKER_LENGTH) {
      bad = 1;
      break;
    }
  }
  if (L > VRNA_GQUAD_MAX_STACK_SIZE)
    bad = 1;

  if (L < VRNA_GQUAD_MIN_STACK_SIZE)
    bad = 1;

  gq = 0.;

  if (!bad) {
    gquad_pf_ali(i, L, l,
                 (void *)(&gq),
                 helper,
                 NULL,
                 NULL);
  }

  for (x = 0; x < L; x++) {
    pp[idx[i + x] - (i + x + 3 * L + l[0] + l[1] + l[2])]                       += gq;
    pp[idx[i + x] - (i + x + L + l[0])]                                         += gq;
    pp[idx[i + x + L + l[0]] - (i + x + 2 * L + l[0] + l[1])]                   += gq;
    pp[idx[i + x + 2 * L + l[0] + l[1]] - (i + x + 3 * L + l[0] + l[1] + l[2])] += gq;
  }
}


PRIVATE INLINE int *
get_g_islands(short *S)
{
  return get_g_islands_sub(S, 1, S[0]);
}


PRIVATE INLINE int *
get_g_islands_sub(short *S,
                  int   i,
                  int   j)
{
  int x, *gg;

  gg  = (int *)vrna_alloc(sizeof(int) * (j - i + 2));
  gg  -= i - 1;

  if (S[j] == 3)
    gg[j] = 1;

  for (x = j - 1; x >= i; x--)
    if (S[x] == 3)
      gg[x] = gg[x + 1] + 1;

  return gg;
}


/**
 *  We could've also created a macro that loops over all G-quadruplexes
 *  delimited by i and j. However, for the fun of it we use this function
 *  that receives a pointer to a callback function which in turn does the
 *  actual computation for each quadruplex found.
 */
PRIVATE
void
process_gquad_enumeration(int *gg,
                          int i,
                          int j,
                          void (*f)(int, int, int *,
                                    void *, void *, void *, void *),
                          void *data,
                          void *P,
                          void *aux1,
                          void *aux2)
{
  int L, l[3], n, max_linker, maxl0, maxl1;

  n = j - i + 1;

  if ((n >= VRNA_GQUAD_MIN_BOX_SIZE) && (n <= VRNA_GQUAD_MAX_BOX_SIZE)) {
    for (L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
         L >= VRNA_GQUAD_MIN_STACK_SIZE;
         L--)
      if (gg[j - L + 1] >= L) {
        max_linker = n - 4 * L;
        if ((max_linker >= 3 * VRNA_GQUAD_MIN_LINKER_LENGTH)
            && (max_linker <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH)) {
          maxl0 = MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH,
                       max_linker - 2 * VRNA_GQUAD_MIN_LINKER_LENGTH
                       );
          for (l[0] = VRNA_GQUAD_MIN_LINKER_LENGTH;
               l[0] <= maxl0;
               l[0]++)
            if (gg[i + L + l[0]] >= L) {
              maxl1 = MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH,
                           max_linker - l[0] - VRNA_GQUAD_MIN_LINKER_LENGTH
                           );
              for (l[1] = VRNA_GQUAD_MIN_LINKER_LENGTH;
                   l[1] <= maxl1;
                   l[1]++)
                if (gg[i + 2 * L + l[0] + l[1]] >= L) {
                  l[2] = max_linker - l[0] - l[1];
                  f(i, L, &(l[0]), data, P, aux1, aux2);
                }
            }
        }
      }
  }
}
