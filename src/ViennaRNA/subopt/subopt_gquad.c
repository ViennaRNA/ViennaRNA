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
#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/subopt/gquad.h"

#include "ViennaRNA/intern/gquad_helpers.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE
void
gquad_pos_exhaustive(unsigned int  i,
                     unsigned int  L,
                     unsigned int  *l,
                     void *data,
                     void *P,
                     void *Lex,
                     void *lex);


PRIVATE
void
gquad_count(unsigned int   i,
            unsigned int   L,
            unsigned int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3);


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */

vrna_array(int)
vrna_gq_int_loop_subopt(vrna_fold_compound_t      *fc,
                        unsigned int              i,
                        unsigned int              j,
                        vrna_array(unsigned int)  *ps,
                        vrna_array(unsigned int)  *qs,
                        int                       threshold){
  short         *S, *S1, si, sj;
  unsigned int  type, p, q, l1, minq, maxq;
  int           energy, *ge, e_gq, dangles, c0;

  ge  = NULL;
  *ps = NULL;
  *qs = NULL;

  if ((fc) &&
      (ps) &&
      (qs)) {
    vrna_param_t  *P  = fc->params;
    vrna_md_t     *md = &(P->model_details);

    vrna_smx_csr(int) * c_gq = fc->matrices->c_gq;

    S       = fc->sequence_encoding2;
    S1      = fc->sequence_encoding;
    type    = vrna_get_ptype_md(S[i], S[j], md);
    dangles = md->dangles;
    si      = S1[i + 1];
    sj      = S1[j - 1];
    energy  = 0;

    if (dangles == 2)
      energy += P->mismatchI[type][si][sj];

    if (type > 2)
      energy += P->TerminalAU;

    vrna_array_init(*ps);
    vrna_array_init(*qs);
    vrna_array_init(ge);

    p = i + 1;
    if (S[p] == 3) {
      if (p + VRNA_GQUAD_MIN_BOX_SIZE < j) {
        minq  = j - i + p - MAXLOOP - 2;
        c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        minq  = MAX2(c0, minq);
        c0    = j - 3;
        maxq  = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        maxq  = MIN2(c0, maxq);
        for (q = minq; q < maxq; q++) {
          if (S[q] != 3)
            continue;

#ifndef VRNA_DISABLE_C11_FEATURES
          e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
          e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

          if (e_gq != INF) {
            c0 = energy + e_gq + P->internal_loop[j - q - 1];

            if (c0 <= threshold) {
              vrna_array_append(ge, energy + P->internal_loop[j - q - 1]);
              vrna_array_append(*ps, p);
              vrna_array_append(*qs, q);
            }
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

#ifndef VRNA_DISABLE_C11_FEATURES
        e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
        e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

        if (e_gq != INF) {
          c0 = energy + e_gq + P->internal_loop[l1 + j - q - 1];

          if (c0 <= threshold) {
            vrna_array_append(ge, energy + P->internal_loop[l1 + j - q - 1]);
            vrna_array_append(*ps, p);
            vrna_array_append(*qs, q);
          }
        }
      }
    }

    q = j - 1;
    if (S[q] == 3)
      for (p = i + 4;
           p + VRNA_GQUAD_MIN_BOX_SIZE < j;
           p++) {
        l1 = p - i - 1;
        if (l1 > MAXLOOP)
          break;

        if (S[p] != 3)
          continue;

#ifndef VRNA_DISABLE_C11_FEATURES
        e_gq = vrna_smx_csr_get(c_gq, p, q, INF);
#else
        e_gq = vrna_smx_csr_int_get(c_gq, p, q, INF);
#endif

        if (e_gq != INF) {
          c0 = energy + e_gq + P->internal_loop[l1];
          if (c0 <= threshold) {
            vrna_array_append(ge, energy + P->internal_loop[l1]);
            vrna_array_append(*ps, p);
            vrna_array_append(*qs, q);
          }
        }
      }
  }

  return ge;
}


PUBLIC void
get_gquad_pattern_exhaustive(short        *S,
                             unsigned int i,
                             unsigned int j,
                             vrna_param_t *P,
                             unsigned int *L,
                             unsigned int *l,
                             int          threshold)
{
  unsigned int *gg;

  gg = get_g_islands_sub(S, i, j);

  process_gquad_enumeration(gg, i, j,
                            &gquad_pos_exhaustive,
                            (void *)(&threshold),
                            (void *)P,
                            (void *)&L,
                            (void *)l);

  gg += i - 1;

  free(gg);
}


PUBLIC unsigned int
get_gquad_count(short         *S,
                unsigned int  i,
                unsigned int  j)
{
  unsigned int p, q, *gg, counter;

  gg = get_g_islands_sub(S, i, j);
  counter = 0;

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


/*
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */
PRIVATE
void
gquad_pos_exhaustive(unsigned int  i,
                     unsigned int  L,
                     unsigned int  *l,
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
    for (cnt = 0; ((unsigned int *)Lex)[cnt] != 0; cnt++);

    *((unsigned int *)Lex + cnt)             = L;
    *((unsigned int *)Lex + cnt + 1)         = 0;
    *(((unsigned int *)lex) + (3 * cnt) + 0) = l[0];
    *(((unsigned int *)lex) + (3 * cnt) + 1) = l[1];
    *(((unsigned int *)lex) + (3 * cnt) + 2) = l[2];
  }
}


PRIVATE
void
gquad_count(unsigned int   i,
            unsigned int   L,
            unsigned int   *l,
            void  *data,
            void  *NA,
            void  *NA2,
            void  *NA3)
{
  *((unsigned int *)data) += 1;
}


