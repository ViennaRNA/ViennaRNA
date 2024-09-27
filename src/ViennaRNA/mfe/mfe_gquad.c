/*
 * G-Quadruplex MFE functions
 *
 * Ronny Lorenz 2024
 *
 * ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/mfe/gquad.h"
#include "ViennaRNA/intern/gquad_helpers.h"


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */


/**
 *  MFE callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe(unsigned int  i,
          unsigned int  L,
          unsigned int  *l,
          void          *data,
          void          *P,
          void          *NA,
          void          *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 */
PRIVATE
void
gquad_mfe_ali(unsigned int  i,
              unsigned int  L,
              unsigned int  *l,
              void          *data,
              void          *helper,
              void          *NA,
              void          *NA2);


/**
 * MFE (alifold) callback for process_gquad_enumeration()
 * with seperation of free energy and penalty contribution
 */
PRIVATE
void
gquad_mfe_ali_en(unsigned int i,
                 unsigned int L,
                 unsigned int *l,
                 void         *data,
                 void         *helper,
                 void         *NA,
                 void         *NA2);


/* other useful static functions */

PRIVATE int
E_gquad_consensus(unsigned int        L,
                  unsigned int        l[3],
                  unsigned int        position,
                  unsigned int        length,
                  unsigned int        n_seq,
                  const unsigned int  **a2s,
                  vrna_param_t        *P);


PRIVATE
int
E_gquad_ali_penalty(unsigned int  L,
                    unsigned int  l[3],
                    unsigned int  i,
                    unsigned int  length,
                    unsigned int  n_seq,
                    const short   **S,
                    vrna_param_t  *P);


PRIVATE void
count_gquad_layer_mismatches(unsigned int L,
                             unsigned int l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2]);


PRIVATE INLINE void
aln_linker_positions(unsigned int L,
                     unsigned int l[3],
                     unsigned int position,
                     unsigned int length,
                     unsigned int starts[3],
                     unsigned int ends[3]);


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


/* retrieve a set of sequence coordinates for the Gs involved
 * in a layer (1-based) of a GQ with stack size L and linker
 * lengths l starting at position i. The GQ may cross the n,1
 * junction so the total length of the sequence (alignment) has
 * to be passed through variable n
 */
PRIVATE void
gq_layer_pos(unsigned int L,
             unsigned int l[3],
             unsigned int layer,
             unsigned int i,
             unsigned int n,
             unsigned int layer_pos[4]);


PRIVATE INLINE int
aln_linker_length(unsigned int        start,
                  unsigned int        end,
                  unsigned int        n,
                  const unsigned int  *a2ss);


PRIVATE void
count_gquad_layer_mismatches(unsigned int L,
                             unsigned int l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2]);


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */
PUBLIC
vrna_smx_csr(int) *
vrna_mfe_gquad_mx(vrna_fold_compound_t * fc){
  vrna_smx_csr(int) * gq_mfe_pos = NULL;

  if (fc) {
    unsigned int  i, j, n, n2;
    unsigned int  *gg;
    vrna_param_t  *P;
    short         *S_enc, *S_tmp;
    void          *data;
    void          ( *process_f )(unsigned int,
                                 unsigned int,
                                 unsigned int *,
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


PUBLIC int
vrna_mfe_gquad_internal_loop(vrna_fold_compound_t *fc,
                             unsigned int         i,
                             unsigned int         j)
{
  unsigned int  type, s, n_seq, **a2s, p, q, u1, u2, minq, maxq, l1, sliding_window;
  int           energy, ge, e_gq, dangles, c0, **ggg;
  short         *S, *S1, si, sj, **SS, **S5, **S3;
  vrna_param_t  *P;
  vrna_md_t     *md;

  vrna_smx_csr(int) * c_gq;

  ge = INF;

  if ((fc) &&
      (i > 0) &&
      (i + VRNA_GQUAD_MIN_BOX_SIZE < j)) {
    sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
    n_seq           = fc->n_seq;
    ggg             = (sliding_window) ? fc->matrices->ggg_local : NULL;
    S               = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : NULL;
    S1              = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding : fc->S_cons;
    SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
    S5              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S5;
    S3              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S3;
    a2s             = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->a2s;
    c_gq            = (sliding_window) ? NULL : fc->matrices->c_gq;
    P               = fc->params;
    md              = &(P->model_details);
    dangles         = md->dangles;
    si              = S1[i + 1];
    sj              = S1[j - 1];
    energy          = 0;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type = vrna_get_ptype_md(S[i], S[j], md);
        if (dangles)
          energy += P->mismatchI[type][si][sj];

        if (type > 2)
          energy += P->TerminalAU;

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        for (s = 0; s < n_seq; s++) {
          type = vrna_get_ptype_md(SS[s][i], SS[s][j], md);
          if (md->dangles)
            energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

          if (type > 2)
            energy += P->TerminalAU;
        }
        break;

      default:
        return ge;
    }


    p = i + 1;
    if (S1[p] == 3) {
      if (p + VRNA_GQUAD_MIN_BOX_SIZE < j) {
        minq = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        if (minq + 1 + MAXLOOP < j)
          minq = j - MAXLOOP - 1;

        maxq = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        if (maxq + 3 > j)
          maxq = j - 3;

        for (q = minq; q < maxq; q++) {
          if (S1[q] != 3)
            continue;

          if (sliding_window) {
            e_gq = ggg[p][q - p];
          } else {
#ifndef VRNA_DISABLE_C11_FEATURES
            e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
            e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif
          }

          if (e_gq != INF) {
            c0 = energy + e_gq;

            switch (fc->type) {
              case VRNA_FC_TYPE_SINGLE:
                c0 += P->internal_loop[j - q - 1];
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                for (s = 0; s < n_seq; s++) {
                  u1  = a2s[s][j - 1] - a2s[s][q];
                  c0  += P->internal_loop[u1];
                }

                break;
            }

            ge = MIN2(ge, c0);
          }
        }
      }
    }

    for (p = i + 2;
         p + VRNA_GQUAD_MIN_BOX_SIZE < j;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S1[p] != 3)
        continue;

      minq = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      if (minq + 1 + MAXLOOP - l1 < j)
        minq = j - MAXLOOP + l1 - 1;

      maxq = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      if (maxq >= j)
        maxq = j - 1;

      for (q = minq; q < maxq; q++) {
        if (S1[q] != 3)
          continue;

        if (sliding_window) {
          e_gq = ggg[p][q - p];
        } else {
#ifndef VRNA_DISABLE_C11_FEATURES
          e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
          e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif
        }

        if (e_gq != INF) {
          c0 = energy + e_gq;

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              c0 += P->internal_loop[l1 + j - q - 1];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                u1  = a2s[s][p - 1] - a2s[s][i];
                u2  = a2s[s][j - 1] - a2s[s][q];
                c0  += P->internal_loop[u1 + u2];
              }
              break;
          }

          ge = MIN2(ge, c0);
        }
      }
    }

    q = j - 1;
    if (S1[q] == 3)
      for (p = (i + 4 + VRNA_GQUAD_MAX_BOX_SIZE - 1 < q) ? q - VRNA_GQUAD_MAX_BOX_SIZE + 1 : i + 4;
           p + VRNA_GQUAD_MIN_BOX_SIZE - 1 < j;
           p++) {
        l1 = p - i - 1;
        if (l1 > MAXLOOP)
          break;

        if (S1[p] != 3)
          continue;

        if (sliding_window) {
          e_gq = ggg[p][q - p];
        } else {
#ifndef VRNA_DISABLE_C11_FEATURES
          e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
          e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif
        }

        if (e_gq != INF) {
          c0 = energy + e_gq;

          switch (fc->type) {
            case VRNA_FC_TYPE_SINGLE:
              c0 += P->internal_loop[l1];
              break;

            case VRNA_FC_TYPE_COMPARATIVE:
              for (s = 0; s < n_seq; s++) {
                u1  = a2s[s][p - 1] - a2s[s][i];
                c0  += P->internal_loop[u1];
              }
              break;
          }

          ge = MIN2(ge, c0);
        }
      }
  }

  return ge;
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
  int           **data;
  unsigned int  i, j, k, *gg, p, q;

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
    for (i = 0; i < (unsigned int)maxdist + 5; i++)
      data[start][i] = INF;

    /*  now we compute contributions for all gquads with 5' delimiter at
     *  position 'start'
     */
    FOR_EACH_GQUAD_AT((unsigned int)start, j, (unsigned int)(start + maxdist + 4)){
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
    for (k = n; ((int)k + maxdist + 5 > n) && (k > 0); k--) {
      data[k] = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
      for (i = 0; i < (unsigned int)maxdist + 5; i++)
        data[k][i] = INF;
    }

    /* compute all contributions for the gquads in this interval */
    unsigned int start = 1;
    if (maxdist + 4 < n)
      start = (unsigned int)(n - maxdist - 4);

    FOR_EACH_GQUAD(i, j, start, (unsigned int)n){
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
  int           **data;
  unsigned int  i, j, k, *gg, p, q;

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
    for (i = 0; i < (unsigned int)maxdist + 5; i++)
      data[start][i] = INF;

    /*  now we compute contributions for all gquads with 5' delimiter at
     *  position 'start'
     */
    FOR_EACH_GQUAD_AT((unsigned int)start, j, (unsigned int)(start + maxdist + 4)){
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
    k = 0;
    if (maxdist + 4 < n)
      k = n - maxdist - 4;

    for (; k <= (unsigned int)n; k++) {
      data[k] = (int *)vrna_alloc(sizeof(int) * (maxdist + 5));
      for (i = 0; i < (unsigned int)maxdist + 5; i++)
        data[k][i] = INF;
    }

    /* compute all contributions for the gquads in this interval */
    unsigned int start = 1;
    if (maxdist + 4 < n)
      start = (unsigned int)(n - maxdist - 4);

    FOR_EACH_GQUAD(i, j, start, (unsigned int)n){
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
PRIVATE void
gquad_mfe(unsigned int  i VRNA_UNUSED,
          unsigned int  L,
          unsigned int  *l,
          void          *data,
          void          *P,
          void          *NA VRNA_UNUSED,
          void          *NA2 VRNA_UNUSED)
{
  int cc = ((vrna_param_t *)P)->gquad[L][l[0] + l[1] + l[2]];

  if (cc < *((int *)data))
    *((int *)data) = cc;
}


PRIVATE void
gquad_mfe_ali(unsigned int  i,
              unsigned int  L,
              unsigned int  *l,
              void          *data,
              void          *helper,
              void          *NA VRNA_UNUSED,
              void          *NA2 VRNA_UNUSED)
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
gquad_mfe_ali_en(unsigned int i,
                 unsigned int L,
                 unsigned int *l,
                 void         *data,
                 void         *helper,
                 void         *NA VRNA_UNUSED,
                 void         *NA2 VRNA_UNUSED)
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


PRIVATE int
E_gquad_consensus(unsigned int        L,
                  unsigned int        l[3],
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


PRIVATE int
E_gquad_ali_penalty(unsigned int  L,
                    unsigned int  l[3],
                    unsigned int  i,
                    unsigned int  length,
                    unsigned int  n_seq,
                    const short   **S,
                    vrna_param_t  *P)
{
  unsigned int mm[2];

  count_gquad_layer_mismatches(L, l, i, length, n_seq, S, mm);

  if (mm[1] > (unsigned int)P->gquadLayerMismatchMax)
    return INF;
  else
    return P->gquadLayerMismatch * (int)mm[0];
}


/* compute (individual) lengths of the unpaired linker sequences */
/*  note here, that we might have a GQ spanning the n,1 junction,
 *  so we first need to transform the linker start and end
 *  positions accordingly
 */
PRIVATE INLINE void
aln_linker_positions(unsigned int L,
                     unsigned int l[3],
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
gq_layer_pos(unsigned int L,
             unsigned int l[3],
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


PRIVATE void
count_gquad_layer_mismatches(unsigned int L,
                             unsigned int l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2])
{
  unsigned int s, layer_pos[4], mismatch, cnt;

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
  unsigned int  n, size, i, j, *gg;

  int           *my_index, *data;

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


PUBLIC int *
get_gquad_ali_matrix(unsigned int n,
                     short        *S_cons,
                     short        **S,
                     unsigned int **a2s,
                     int          n_seq,
                     vrna_param_t *P)
{
  unsigned int  i, j, size, *gg;
  int           *my_index, *data;

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


PUBLIC int
E_GQuad_IntLoop_L_comparative(int           i,
                              int           j,
                              unsigned int  *tt,
                              short         *S_cons,
                              short         **S5,
                              short         **S3,
                              unsigned int  **a2s,
                              int           **ggg,
                              int           n_seq,
                              vrna_param_t  *P)
{
  unsigned int  type;
  int           eee, energy, ge, p, q, l1, u1, u2, minq, maxq, c0, s;
  vrna_md_t     *md;

  md      = &(P->model_details);
  energy  = 0;

  for (s = 0; s < n_seq; s++) {
    type = tt[s];
    if (md->dangles == 2)
      energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

    if (type > 2)
      energy += P->TerminalAU;
  }

  ge = INF;

  p = i + 1;
  if (S_cons[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S_cons[q] != 3)
          continue;

        eee = 0;

        for (s = 0; s < n_seq; s++) {
          u1  = a2s[s][j - 1] - a2s[s][q];
          eee += P->internal_loop[u1];
        }

        c0 = energy +
             ggg[p][q - p] +
             eee;
        ge = MIN2(ge, c0);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S_cons[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S_cons[q] != 3)
        continue;

      eee = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][p - 1] - a2s[s][i];
        u2  = a2s[s][j - 1] - a2s[s][q];
        eee += P->internal_loop[u1 + u2];
      }

      c0 = energy +
           ggg[p][q - p] +
           eee;
      ge = MIN2(ge, c0);
    }
  }

  q = j - 1;
  if (S_cons[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S_cons[p] != 3)
        continue;

      eee = 0;

      for (s = 0; s < n_seq; s++) {
        u1  = a2s[s][p - 1] - a2s[s][i];
        eee += P->internal_loop[u1];
      }

      c0 = energy +
           ggg[p][q - p] +
           eee;
      ge = MIN2(ge, c0);
    }

  return ge;
}


PUBLIC int
E_GQuad_IntLoop_L(int           i,
                  int           j,
                  int           type,
                  short         *S,
                  int           **ggg,
                  int           maxdist VRNA_UNUSED,
                  vrna_param_t  *P)
{
  int   energy, ge, dangles, p, q, l1, minq, maxq, c0;
  short si, sj;

  dangles = P->model_details.dangles;
  si      = S[i + 1];
  sj      = S[j - 1];
  energy  = 0;

  if (dangles == 2)
    energy += P->mismatchI[type][si][sj];

  if (type > 2)
    energy += P->TerminalAU;

  ge = INF;

  p = i + 1;
  if (S[p] == 3) {
    if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
      minq  = j - i + p - MAXLOOP - 2;
      c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
      minq  = MAX2(c0, minq);
      c0    = j - 3;
      maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
      maxq  = MIN2(c0, maxq);
      for (q = minq; q < maxq; q++) {
        if (S[q] != 3)
          continue;

        c0  = energy + ggg[p][q - p] + P->internal_loop[j - q - 1];
        ge  = MIN2(ge, c0);
      }
    }
  }

  for (p = i + 2;
       p < j - VRNA_GQUAD_MIN_BOX_SIZE;
       p++) {
    l1 = p - i - 1;
    if (l1 > MAXLOOP)
      break;

    if (S[p] != 3)
      continue;

    minq  = j - i + p - MAXLOOP - 2;
    c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
    minq  = MAX2(c0, minq);
    c0    = j - 1;
    maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
    maxq  = MIN2(c0, maxq);
    for (q = minq; q < maxq; q++) {
      if (S[q] != 3)
        continue;

      c0  = energy + ggg[p][q - p] + P->internal_loop[l1 + j - q - 1];
      ge  = MIN2(ge, c0);
    }
  }

  q = j - 1;
  if (S[q] == 3)
    for (p = i + 4;
         p < j - VRNA_GQUAD_MIN_BOX_SIZE;
         p++) {
      l1 = p - i - 1;
      if (l1 > MAXLOOP)
        break;

      if (S[p] != 3)
        continue;

      c0  = energy + ggg[p][q - p] + P->internal_loop[l1];
      ge  = MIN2(ge, c0);
    }

  return ge;
}


#endif
