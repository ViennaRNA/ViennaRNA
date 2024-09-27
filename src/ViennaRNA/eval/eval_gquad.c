/*
 * eval/eval_gquad.c
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
#include <math.h>
#include <string.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/eval/gquad.h"

#include "ViennaRNA/intern/gquad_helpers.h"


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

PRIVATE
int
E_gquad_ali_penalty(unsigned int  L,
                    unsigned int  l[3],
                    unsigned int  i,
                    unsigned int  length,
                    unsigned int  n_seq,
                    const short   **S,
                    vrna_param_t  *P);


PRIVATE
FLT_OR_DBL
exp_E_gquad_ali_penalty(unsigned int      L,
                        unsigned int      l[3],
                        unsigned int      i,
                        unsigned int      length,
                        unsigned int      n_seq,
                        const short       **S,
                        vrna_exp_param_t  *P);


PRIVATE int
E_gquad_consensus(unsigned int        L,
                  unsigned int        l[3],
                  unsigned int        position,
                  unsigned int        length,
                  unsigned int        n_seq,
                  const unsigned int  **a2s,
                  vrna_param_t        *P);


PRIVATE FLT_OR_DBL
exp_E_gquad_consensus(unsigned int        L,
                      unsigned int        l[3],
                      unsigned int        position,
                      unsigned int        length,
                      unsigned int        n_seq,
                      const unsigned int  **a2s,
                      vrna_exp_param_t    *pf);


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


PRIVATE INLINE int
aln_linker_length(unsigned int        start,
                  unsigned int        end,
                  unsigned int        n,
                  const unsigned int  *a2ss);


PRIVATE void
gq_layer_pos(unsigned int L,
             unsigned int l[3],
             unsigned int layer,
             unsigned int i,
             unsigned int n,
             unsigned int layer_pos[4]);


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */
PUBLIC int
vrna_E_gquad(unsigned int L,
             unsigned int l[3],
             vrna_param_t *P)
{
  if (P) {
    CHECK_GQUAD(L, l, return INF);

    return P->gquad[L][l[0] + l[1] + l[2]];
  }

  return INF;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_gquad(unsigned int     L,
                 unsigned int     l[3],
                 vrna_exp_param_t *pf)
{
  if (pf) {
    CHECK_GQUAD(L, l, return 0.);

    return pf->expgquad[L][l[0] + l[1] + l[2]];
  }

  return 0.;
}


PUBLIC void
vrna_E_consensus_gquad(unsigned int       L,
                       unsigned int       l[3],
                       unsigned int       position,
                       unsigned int       length,
                       unsigned int       n_seq,
                       const short        **S,
                       const unsigned int **a2s,
                       vrna_param_t       *P,
                       int                en[2])
{
  int penalty;

  en[0] = en[1] = INF;

  if (P) {
    CHECK_GQUAD(L, l, return );

    /* get penalty from incompatible sequences in alignment */
    penalty = E_gquad_ali_penalty(L,
                                  l,
                                  position,
                                  length,
                                  n_seq,
                                  S,
                                  P);

    /* assign return values */
    if (penalty != INF) {
      /* compute actual quadruplex contribution for subalignment */
      en[0] = E_gquad_consensus(L,
                                l,
                                position,
                                length,
                                n_seq,
                                (const unsigned int **)a2s,
                                P);
      en[1] = penalty;
    }
  }
}


PUBLIC FLT_OR_DBL
vrna_exp_E_consensus_gquad(unsigned int       L,
                           unsigned int       l[3],
                           vrna_exp_param_t   *pf,
                           unsigned int       position,
                           unsigned int       length,
                           unsigned int       n_seq,
                           const short        **S,
                           const unsigned int **a2s)
{
  FLT_OR_DBL penalty, q = 0.;

  if ((pf) &&
      (S) &&
      (a2s)) {
    CHECK_GQUAD(L, l, return q);

    penalty = exp_E_gquad_ali_penalty(L,
                                      l,
                                      position,
                                      length,
                                      n_seq,
                                      S,
                                      pf);

    if (penalty != 0.)
      q = penalty *
          exp_E_gquad_consensus(L,
                                l,
                                position,
                                length,
                                n_seq,
                                a2s,
                                pf);
  }

  return q;
}


/*
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */
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

  if (mm[1] > P->gquadLayerMismatchMax)
    return INF;
  else
    return P->gquadLayerMismatch * mm[0];
}


PRIVATE FLT_OR_DBL
exp_E_gquad_ali_penalty(unsigned int      L,
                        unsigned int      l[3],
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


PRIVATE FLT_OR_DBL
exp_E_gquad_consensus(unsigned int        L,
                      unsigned int        l[3],
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


PRIVATE void
count_gquad_layer_mismatches(unsigned int L,
                             unsigned int l[3],
                             unsigned int i,
                             unsigned int n,
                             unsigned int n_seq,
                             const short  **S,
                             unsigned int mm[2])
{
  unsigned int  s, layer_pos[4], mismatch, cnt;

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


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PUBLIC int
E_gquad(int           L,
        int           l[3],
        vrna_param_t  *P)
{
  unsigned int LL, ll[3];

  LL    = (unsigned int)L;
  ll[0] = (unsigned int)l[0];
  ll[1] = (unsigned int)l[1];
  ll[2] = (unsigned int)l[2];

  return vrna_E_gquad(LL, ll, P);
}


PUBLIC FLT_OR_DBL
exp_E_gquad(int               L,
            int               l[3],
            vrna_exp_param_t  *pf)
{
  unsigned int LL, ll[3];

  LL    = (unsigned int)L;
  ll[0] = (unsigned int)l[0];
  ll[1] = (unsigned int)l[1];
  ll[2] = (unsigned int)l[2];

  return vrna_exp_E_gquad(LL, ll, pf);
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
  unsigned int LL, ll[3];

  LL    = (unsigned int)L;
  ll[0] = (unsigned int)l[0];
  ll[1] = (unsigned int)l[1];
  ll[2] = (unsigned int)l[2];

  vrna_E_consensus_gquad(LL,
                         ll,
                         (unsigned int)i,
                         0,
                         n_seq,
                         S,
                         (const unsigned int **)a2s,
                         P,
                         en);
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
  unsigned int LL, ll[3];

  LL    = (unsigned int)L;
  ll[0] = (unsigned int)l[0];
  ll[1] = (unsigned int)l[1];
  ll[2] = (unsigned int)l[2];

  return vrna_exp_E_consensus_gquad(LL,
                                    ll,
                                    pf,
                                    (unsigned int)i,
                                    0,
                                    n_seq,
                                    (const short **)S,
                                    (const unsigned int **)a2s);
}


#endif
