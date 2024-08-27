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
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/loops/gquad.h"

#include "ViennaRNA/loops/gquad_intern.h"


#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

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
                          void ( *f )(int, int, int *,
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


/* other useful static functions */

PRIVATE int
E_gquad_consensus(int                 L,
                  int                 l[3],
                  unsigned int        position,
                  unsigned int        length,
                  unsigned int        n_seq,
                  const unsigned int  **a2s,
                  vrna_param_t        *P);


PRIVATE
int
E_gquad_ali_penalty(int           L,
                    int           l[3],
                    unsigned int  i,
                    unsigned int  length,
                    unsigned int  n_seq,
                    const short   **S,
                    vrna_param_t  *P);


PRIVATE
FLT_OR_DBL
exp_E_gquad_ali_penalty(int               L,
                        int               l[3],
                        unsigned int      i,
                        unsigned int      length,
                        unsigned int      n_seq,
                        const short       **S,
                        vrna_exp_param_t  *P);


PRIVATE void
count_gquad_layer_mismatches(int          L,
                             int          l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
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


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */

/********************************
 * Here are the G-quadruplex
 * dynamic programming matrices
 ********************************/

PUBLIC
vrna_smx_csr(int) *
vrna_gq_pos_mfe(vrna_fold_compound_t * fc){
  vrna_smx_csr(int) * gq_mfe_pos = NULL;

  if (fc) {
    int           i, j, n, n2;
    int           *gg;
    vrna_param_t  *P;
    short         *S_enc, *S_tmp;
    void          *data;
    void          ( *process_f )(int,
                                 int,
                                 int *,
                                 void *,
                                 void *,
                                 void *,
                                 void *);
    struct gquad_ali_helper tmp = {
      0
    };

    n           = fc->length;
    n2          = 0;
    P           = fc->params;
    S_tmp       = NULL;
    gq_mfe_pos  = vrna_smx_csr_int_init(n + 1);

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S_enc     = fc->sequence_encoding2;
        data      = (void *)P;
        process_f = &gquad_mfe;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        S_enc = fc->S_cons;
        struct gquad_ali_helper gq_help = {
          .S      = (const short **)fc->S,
          .a2s    = (const unsigned int **)fc->a2s,
          .length = fc->length,
          .n_seq  = fc->n_seq,
          .P      = P
        };
        tmp       = gq_help;
        data      = (void *)&tmp;
        process_f = &gquad_mfe_ali;
        break;

      default:
        return NULL;
    }

    if (P->model_details.circ) {
      n2 = MIN2(n, VRNA_GQUAD_MAX_BOX_SIZE) - 1;

      S_tmp = (short *)vrna_alloc(sizeof(short) * (n + n2 + 1));
      memcpy(S_tmp, S_enc, sizeof(short) * (n + 1));
      memcpy(S_tmp + (n + 1), S_enc + 1, sizeof(short) * n2);
      S_tmp[0]  = n + n2;
      S_enc     = S_tmp;
      n         += n2;
    }

    gg = get_g_islands(S_enc);

    FOR_EACH_GQUAD_INC(i, j, 1, n) {
      int e = INF;

      if (i > n - n2)
        break;

      process_gquad_enumeration(gg, i, j,
                                process_f,
                                (void *)(&e),
                                data,
                                NULL,
                                NULL);
      if ((e < INF) && (j - i + 1 <= n - n2)) {
#ifndef VRNA_DISABLE_C11_FEATURES
        vrna_smx_csr_insert(gq_mfe_pos, i, (j - 1) % (n - n2) + 1, e);
#else
        vrna_smx_csr_int_insert(gq_mfe_pos, i, (j - 1) % (n - n2) + 1, e);
#endif
      }
    }

    free(S_tmp);
    free(gg);
  }

  return gq_mfe_pos;
}


PUBLIC
vrna_smx_csr(FLT_OR_DBL) *
vrna_gq_pos_pf(vrna_fold_compound_t * fc){
  vrna_smx_csr(FLT_OR_DBL) * q_gq = NULL;

  if (fc) {
    int               i, j, n, n2, *gg;
    short             *S_tmp, *S_enc;
    FLT_OR_DBL        q, *scale;
    vrna_exp_param_t  *pf_params;
    void              *data;
    void              ( *process_f )(int,
                                     int,
                                     int *,
                                     void *,
                                     void *,
                                     void *,
                                     void *);
    struct gquad_ali_helper tmp = {
      0
    };

    n         = fc->length;
    n2        = 0;
    pf_params = fc->exp_params;
    S_tmp     = NULL;
    q_gq      = vrna_smx_csr_FLT_OR_DBL_init(n + 1);
    scale     = fc->exp_matrices->scale;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S_enc     = fc->sequence_encoding2;
        data      = (void *)pf_params;
        process_f = &gquad_pf;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        S_enc = fc->S_cons;
        struct gquad_ali_helper gq_help = {
          .S      = (const short **)fc->S,
          .a2s    = (const unsigned int **)fc->a2s,
          .length = fc->length,
          .n_seq  = fc->n_seq,
          .pf     = pf_params
        };
        tmp       = gq_help;
        data      = (void *)&tmp;
        process_f = &gquad_pf_ali;
        break;

      default:
        return NULL;
    }

    if (pf_params->model_details.circ) {
      n2 = MIN2(n, VRNA_GQUAD_MAX_BOX_SIZE) - 1;

      S_tmp = (short *)vrna_alloc(sizeof(short) * (n + n2 + 1));
      memcpy(S_tmp, S_enc, sizeof(short) * (n + 1));
      memcpy(S_tmp + (n + 1), S_enc + 1, sizeof(short) * n2);
      S_tmp[0]  = n + n2;
      S_enc     = S_tmp;
      n         += n2;
    }

    gg = get_g_islands(S_enc);

    FOR_EACH_GQUAD_INC(i, j, 1, n) {
      q = 0.;

      if (i > n - n2)
        break;

      process_gquad_enumeration(gg, i, j,
                                process_f,
                                (void *)(&q),
                                data,
                                NULL,
                                NULL);
      if ((q != 0.) &&
          (j - i + 1 <= n - n2)) {
#ifndef VRNA_DISABLE_C11_FEATURES
        vrna_smx_csr_insert(q_gq, i, (j - 1) % (n - n2) + 1, q * scale[j - i + 1]);

#else
        vrna_smx_csr_FLT_OR_DBL_insert(q_gq, i, (j - 1) % (n - n2) + 1, q * scale[j - i + 1]);
#endif
      }
    }

    free(S_tmp);
    free(gg);
  }

  return q_gq;
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

  struct gquad_ali_helper gq_help = {
    .S      = (const short **)S,
    .a2s    = (const unsigned int **)a2s,
    .length = n,
    .n_seq  = n_seq,
    .P      = P
  };

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


/*
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */
PRIVATE INLINE int
aln_linker_length(unsigned int        start,
                  unsigned int        end,
                  unsigned int        n,
                  const unsigned int  *a2ss)
{
  if (start <= end) {
    return a2ss[end] - a2ss[start - 1];
  } else {
    return a2ss[n] - a2ss[start - 1] + a2ss[end];
  }
}


/* compute (individual) lengths of the unpaired linker sequences */
/*  note here, that we might have a GQ spanning the n,1 junction,
 *  so we first need to transform the linker start and end
 *  positions accordingly
 */
PRIVATE INLINE void
aln_linker_positions(int          L,
                     int          l[3],
                     unsigned int position,
                     unsigned int length,
                     unsigned int starts[3],
                     unsigned int ends[3])
{
  if ((length > 0) &&
      (position + 4 * L + l[0] + l[1] + l[2] >= length)) {
    starts[0] = (position + L - 1) % (length) + 1;
    ends[0]   = (position + L + l[0] - 1 - 1) % (length) + 1;
    starts[1] = (position + 2 * L + l[0] - 1) % (length) + 1;
    ends[1]   = (position + 2 * L + l[0] + l[1] - 1 - 1) % (length) + 1;
    starts[2] = (position + 3 * L + l[0] + l[1] - 1) % (length) + 1;
    ends[2]   = (position + 3 * L + l[0] + l[1] + l[2] - 1 - 1) % (length) + 1;
  } else {
    starts[0] = position + L;
    ends[0]   = starts[0] + l[0] - 1;
    starts[1] = ends[0] + L + 1;
    ends[1]   = starts[1] + l[1] - 1;
    starts[2] = ends[1] + L + 1;
    ends[2]   = starts[2] + l[2] - 1;
  }
}


/* retrieve a set of sequence coordinates for the Gs involved
 * in a layer (1-based) of a GQ with stack size L and linker
 * lengths l starting at position i. The GQ may cross the n,1
 * junction so the total length of the sequence (alignment) has
 * to be passed through variable n
 */
PRIVATE void
gq_layer_pos(int          L,
             int          l[3],
             unsigned int layer,
             unsigned int i,
             unsigned int n,
             unsigned int layer_pos[4])
{
  if ((n > 0) &&
      (i + 4 * L + l[0] + l[1] + l[2] >= n)) {
    layer_pos[0]  = (i + layer - 1 - 1) % (n) + 1;
    layer_pos[1]  = (i + layer + L + l[0] - 1 - 1) % (n) + 1;
    layer_pos[2]  = (i + layer + 2 * L + l[0] + l[1] - 1 - 1) % (n) + 1;
    layer_pos[3]  = (i + layer + 3 * L + l[0] + l[1] + l[2] - 1 - 1) % (n) + 1;
  } else {
    layer_pos[0]  = i + layer - 1;
    layer_pos[1]  = i + layer + L + l[0] - 1;
    layer_pos[2]  = i + layer + 2 * L + l[0] + l[1] - 1;
    layer_pos[3]  = i + layer + 3 * L + l[0] + l[1] + l[2] - 1;
  }
}


PRIVATE int
E_gquad_consensus(int                 L,
                  int                 l[3],
                  unsigned int        position,
                  unsigned int        length,
                  unsigned int        n_seq,
                  const unsigned int  **a2s,
                  vrna_param_t        *P)
{
  unsigned int  l_start[3], l_end[3], s, u1, u2, u3;
  int           e;


  e = 0;

  aln_linker_positions(L, l, position, length, l_start, l_end);

  for (s = 0; s < n_seq; s++) {
    u1  = aln_linker_length(l_start[0], l_end[0], length, a2s[s]);
    u2  = aln_linker_length(l_start[1], l_end[1], length, a2s[s]);
    u3  = aln_linker_length(l_start[2], l_end[2], length, a2s[s]);

    e += P->gquad[L][u1 + u2 + u3];
  }

  return e;
}


PRIVATE FLT_OR_DBL
exp_E_gquad_consensus(int                 L,
                      int                 l[3],
                      unsigned int        position,
                      unsigned int        length,
                      unsigned int        n_seq,
                      const unsigned int  **a2s,
                      vrna_exp_param_t    *pf)
{
  unsigned int  l_start[3], l_end[3], s, u1, u2, u3;
  FLT_OR_DBL    q;


  q = 1.;

  aln_linker_positions(L, l, position, length, l_start, l_end);

  for (s = 0; s < n_seq; s++) {
    u1  = aln_linker_length(l_start[0], l_end[0], length, a2s[s]);
    u2  = aln_linker_length(l_start[1], l_end[1], length, a2s[s]);
    u3  = aln_linker_length(l_start[2], l_end[2], length, a2s[s]);

    q *= pf->expgquad[L][u1 + u2 + u3];
  }

  return q;
}


PRIVATE int
E_gquad_ali_penalty(int           L,
                    int           l[3],
                    unsigned int  i,
                    unsigned int  length,
                    unsigned int  n_seq,
                    const short   **S,
                    vrna_param_t  *P)
{
  unsigned int mm[2];

  count_gquad_layer_mismatches(L, l, i, length, n_seq, S, mm);

  if (mm[1] > P->gquadLayerMismatchMax)
    return INF;
  else
    return P->gquadLayerMismatch * mm[0];
}


PRIVATE FLT_OR_DBL
exp_E_gquad_ali_penalty(int               L,
                        int               l[3],
                        unsigned int      i,
                        unsigned int      n,
                        unsigned int      n_seq,
                        const short       **S,
                        vrna_exp_param_t  *pf)
{
  unsigned int mm[2];

  count_gquad_layer_mismatches(L, l, i, n, n_seq, S, mm);

  if (mm[1] > pf->gquadLayerMismatchMax)
    return (FLT_OR_DBL)0.;
  else
    return (FLT_OR_DBL)pow(pf->expgquadLayerMismatch, (double)mm[0]);
}


PRIVATE void
count_gquad_layer_mismatches(int          L,
                             int          l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2])
{
  unsigned int  s, layer_pos[4], mismatch;
  int           cnt;

  mm[0] = mm[1] = 0;


  /* check for compatibility in the alignment */
  for (s = 0; s < n_seq; s++) {
    mismatch = 0;

    /* check bottom layer */
    gq_layer_pos(L, l, 1, i, n, layer_pos);

    if (S[s][layer_pos[0]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[1]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[2]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[3]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    }

    /* check top layer */
    gq_layer_pos(L, l, L, i, n, layer_pos);

    if (S[s][layer_pos[0]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[1]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[2]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    } else if (S[s][layer_pos[3]] != 3) {
      /* add 1x penalty for missing bottom layer */
      mismatch++;
    }

    /* check inner layers */
    for (cnt = 2; cnt < L; cnt++) {
      gq_layer_pos(L, l, cnt, i, n, layer_pos);

      if (S[s][layer_pos[0]] != 3) {
        /* add 2x penalty for missing inner layer */
        mismatch += 2;
      } else if (S[s][layer_pos[1]] != 3) {
        /* add 2x penalty for missing inner layer */
        mismatch += 2;
      } else if (S[s][layer_pos[2]] != 3) {
        /* add 2x penalty for missing inner layer */
        mismatch += 2;
      } else if (S[s][layer_pos[3]] != 3) {
        /* add 2x penalty for missing inner layer */
        mismatch += 2;
      }
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
  const short             **S;
  const unsigned int      **a2s;
  unsigned int            n, n_seq;
  FLT_OR_DBL              penalty;
  vrna_exp_param_t        *pf;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n       = gq_help->length;
  n_seq   = gq_help->n_seq;
  pf      = gq_help->pf;
  penalty = exp_E_gquad_ali_penalty(L, l, (unsigned int)i, n, n_seq, S, pf);

  if (penalty != 0.)
    *((FLT_OR_DBL *)data) += penalty *
                             exp_E_gquad_consensus(L, l, (unsigned int)i, n, n_seq, a2s, pf);
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
  int en[2], cc;

  en[0] = en[1] = INF;

  CHECK_GQUAD(L, l, return );

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
  const short             **S;
  const unsigned int      **a2s;
  unsigned int            n_seq, n;
  int                     en[2], cc, dd;
  vrna_param_t            *P;
  struct gquad_ali_helper *gq_help;

  gq_help = (struct gquad_ali_helper *)helper;
  S       = gq_help->S;
  a2s     = gq_help->a2s;
  n       = gq_help->length;
  n_seq   = gq_help->n_seq;
  P       = gq_help->P;

  en[0] = E_gquad_consensus(L, l, i, n, n_seq, (const unsigned int **)a2s, P);
  en[1] = E_gquad_ali_penalty(L, l, i, n, n_seq, S, P);

  if (en[1] != INF) {
    cc  = en[0] + en[1];
    dd  = ((int *)data)[0] + ((int *)data)[1];
    if (cc < dd) {
      ((int *)data)[0]  = en[0];
      ((int *)data)[1]  = en[1];
    }
  }
}


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

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
get_gquad_pf_matrix_comparative(unsigned int      n,
                                short             *S_cons,
                                short             **S,
                                unsigned int      **a2s,
                                FLT_OR_DBL        *scale,
                                unsigned int      n_seq,
                                vrna_exp_param_t  *pf)
{
  int         size, *gg, i, j, *my_index;
  FLT_OR_DBL  *data;


  size      = (n * (n + 1)) / 2 + 2;
  data      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
  gg        = get_g_islands(S_cons);
  my_index  = vrna_idx_row_wise(n);

  struct gquad_ali_helper gq_help = {
    .S      = (const short **)S,
    .a2s    = (const unsigned int **)a2s,
    .length = n,
    .n_seq  = n_seq,
    .pf     = pf
  };

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
get_gquad_ali_matrix(unsigned int n,
                     short        *S_cons,
                     short        **S,
                     unsigned int **a2s,
                     int          n_seq,
                     vrna_param_t *P)
{
  int size, *data, *gg;
  int i, j, *my_index;

  size      = (n * (n + 1)) / 2 + 2;
  data      = (int *)vrna_alloc(sizeof(int) * size);
  gg        = get_g_islands(S_cons);
  my_index  = vrna_idx_col_wise(n);

  struct gquad_ali_helper gq_help = {
    .S      = (const short **)S,
    .a2s    = (const unsigned int **)a2s,
    .length = n,
    .n_seq  = n_seq,
    .P      = P
  };

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


#endif
