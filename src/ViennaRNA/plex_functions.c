/*
 *         compute potentially pseudoknotted structures of an RNA sequence
 *                           Ivo Hofacker
 *                        Vienna RNA package
 */

/*
 * library containing the function used in PKplex
 * it generates pseudoknotted structures by letting the sequence form a duplex structure with itself
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <time.h>

#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/alphabet.h"

#include "ViennaRNA/loops/external_hc.inc"

#include "ViennaRNA/PKplex.h"

#undef  MAXLOOP
#define MAXLOOP 10


PRIVATE vrna_pkplex_t *
duplexfold_XS(vrna_fold_compound_t *fc,
              const int     **access_s1,
              const int     penalty,
              const int     max_interaction_length);


PRIVATE char *
backtrack_XS(vrna_fold_compound_t *fc,
             int        kk,
             int        ll,
             const int  ii,
             const int  jj,
             const int  max_interaction_length,
             int        ***c3);


PRIVATE int
prepare_PKplex(vrna_fold_compound_t *fc);


PUBLIC vrna_pkplex_t *
vrna_PKplex(vrna_fold_compound_t  *fc,
            const int             **accessibility,
            int                   penalty,
            int                   delta,
            unsigned int          max_interaction_length,
            unsigned int          options)
{
  if (fc) {
    prepare_PKplex(fc);

    return duplexfold_XS(fc,
                         accessibility,
                         penalty,
                         max_interaction_length);
  }

  return NULL;
}

PUBLIC vrna_pkplex_t *
PKLduplexfold_XS(const char *s1,
                 const int  **access_s1,
                 int  penalty,
                 int  max_interaction_length,
                 int  delta)
{
  vrna_fold_compound_t  *fc;
  vrna_pkplex_t         *hits;

  fc = vrna_fold_compound(s1, NULL, VRNA_OPTION_DEFAULT);

  prepare_PKplex(fc);

  hits = duplexfold_XS(fc,
                       access_s1,
                       penalty,
                       max_interaction_length);

  vrna_fold_compound_free(fc);

  return hits;
}


PRIVATE int ***
get_array(unsigned int n,
          unsigned int interaction_length)
{
  unsigned int  i, j, k;
  int           ***c3;

  c3 = (int ***)vrna_alloc(sizeof(int **) * (n));

  for (i = 0; i < n; i++) {
    c3[i] = (int **)vrna_alloc(sizeof(int *) * interaction_length);

    for (j = 0; j < interaction_length; j++)
      c3[i][j] = (int *)vrna_alloc(sizeof(int) * interaction_length);
  }

  return c3;
}


PRIVATE void
reset_array(int ***c3,
            unsigned int n,
            unsigned int interaction_length)
{
  unsigned int  i, j, k;

  for (i = 0; i < n; i++) {
    for (j = 0; j < interaction_length; j++) {
      for (k = 0; k < interaction_length; k++)
        c3[i][j][k] = INF;
    }
  }
}


PRIVATE void
free_array(int ***c3,
           unsigned int n,
           unsigned int interaction_length)
{
  unsigned int i, j;

  for (i = 0; i < n; i++) {
    for (j = 0; j < interaction_length; j++)
      free(c3[i][j]);
    free(c3[i]);
  }

  free(c3);
}


PRIVATE int
prepare_PKplex(vrna_fold_compound_t *fc)
{
  /* prepare Boltzmann factors if required */
  vrna_params_prepare(fc, VRNA_OPTION_MFE);

  /* prepare ptype array(s) */
  vrna_ptypes_prepare(fc, VRNA_OPTION_MFE);

  /* prepare hard constraints */
  vrna_hc_prepare(fc, VRNA_OPTION_MFE);

  /* prepare soft constraints data structure, if required */
  vrna_sc_prepare(fc, VRNA_OPTION_MFE);
}


PRIVATE vrna_pkplex_t *
duplexfold_XS(vrna_fold_compound_t *fc,
              const int   **access_s1,
              const int   penalty,
              const int   max_interaction_length)
{
  char          *struc;
  short         *S, *S1, si, sk, sl, sp, sq;
  size_t        storage_size, storage_fill;
  unsigned int  n, type, type2, type3;
  int           ***c3, i, j, k, l, p, q, Emin, l_min, k_min, j_min, E,
                tempK, *rtype, i_pos_begin, j_pos_end, dGx, dGy, inter, turn;

  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_pkplex_t *storage;

  vrna_callback_hc_evaluate *evaluate_ext;
  struct default_data       hc_dat_local;

  struc = NULL;
  n     = fc->length;
  S     = fc->sequence_encoding2;
  S1    = fc->sequence_encoding;
  P     = fc->params;
  md    = &(P->model_details);
  turn  = md->min_loop_size;
  rtype = &(md->rtype[0]);
  hc    = fc->hc;

  storage_size  = 64;
  storage_fill  = 0;
  storage       = (vrna_pkplex_t *)vrna_alloc(sizeof(vrna_pkplex_t) * storage_size);

  evaluate_ext  = prepare_hc_default(fc, &hc_dat_local);

  c3  = get_array(n, max_interaction_length);

  if (n > turn + 1) {
    for (i = n - turn - 1; i > 0; i--) {
      Emin  = INF;
      j_min = 0;
      l_min = 0;
      k_min = 0;

      reset_array(c3, n, max_interaction_length);

      si = S1[i + 1];

      /* matrix starting values for (i,j)-basepairs */
      for (j = i + turn + 1; j <= n; j++) {
        if (evaluate_ext(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          type = md->pair[S[j]][S[i]];
          c3[j - 1][max_interaction_length - 1][0] = vrna_E_ext_stem(type,
                                                                 S1[j - 1],
                                                                 si,
                                                                 P);
        }
      }

      i_pos_begin = MAX2(0, i - max_interaction_length);

      /* fill matrix */
      for (k = i - 1; k > i_pos_begin; k--) {
        tempK = max_interaction_length - i + k - 1;
        sk    = S1[k + 1];
        for (l = i + turn + 1; l <= n; l++) {
          if (hc->mx[n * k + l] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
            /* again, why 9 less then the sequence length ? */
            type2 = md->pair[S[k]][S[l]];
            sl    = S1[l - 1];

            for (p = k + 1; (p <= i) && (p <= k + MAXLOOP + 1); p++) {
              sp  = S1[p - 1];

              for (q = l - 1; (q >= i + turn + 1) && (q >= l - MAXLOOP - 1); q--) {
                if (p - k + l - q - 2 > MAXLOOP)
                  break;

                if (hc->mx[n * p + q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
                  type3 = md->pair[S[q]][S[p]];
                  sq    = S1[q + 1];

                  E = E_IntLoop(p - k - 1,
                                l - q - 1,
                                type2,
                                type3,
                                sk,
                                sl,
                                sp,
                                sq,
                                P);
                  for (j = MAX2(i + turn + 1, l - max_interaction_length + 1); j <= q; j++) {
                    if (hc->mx[n * i + j] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
                      type = md->pair[S[i]][S[j]];
                      c3[j - 1][tempK][l - j] =
                          MIN2(c3[j - 1][tempK][l - j],
                               c3[j - 1][max_interaction_length - i + p - 1][q - j] + E);
                    }
                  }
                }
              } /* next j */
            }   /* next q */
          }     /* next p */
        }       /* next l */
      }         /* next k */

      /* read out matrix minimum */
      for (j = i + turn + 1; j <= n; j++) {
        if (evaluate_ext(i, j, i, j, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
          j_pos_end = MIN2(n + 1, j + max_interaction_length);
          for (k = i - 1; k > i_pos_begin; k--) {
            sk = (k > i_pos_begin + 1) ? S1[k - 1] : -1; /* should actually be k > 10 */
            for (l = j + 1; l < j_pos_end; l++) {
              if (evaluate_ext(k, l, k, l, VRNA_DECOMP_EXT_STEM, &hc_dat_local)) {
                type2 = md->pair[S[k]][S[l]];
                sl  = (l < j_pos_end - 1) ? S1[l + 1] : -1; /* should actually be l < n - 10 */
                E   = c3[j - 1][max_interaction_length - i + k - 1][l - j] +
                      vrna_E_ext_stem(type2, sk, sl, P) +
                      access_s1[i - k + 1][i] +
                      access_s1[l - j + 1][l];

                if (E < Emin) {
                  Emin  = E;
                  k_min = k;
                  l_min = l;
                  j_min = j;
                }
              }
            }
          }
        }
      }

      if (Emin < penalty) {
        struc = backtrack_XS(fc, k_min, l_min, i, j_min, max_interaction_length, c3);

        dGx   = access_s1[i - k_min + 1][i];
        dGy   = access_s1[l_min - j_min + 1][l_min];
        inter = Emin - dGx - dGy;

        storage[storage_fill].tb        = k_min;
        storage[storage_fill].te        = i;
        storage[storage_fill].qb        = j_min;
        storage[storage_fill].qe        = l_min;
        storage[storage_fill].ddG       = (double)Emin * 0.01;
        storage[storage_fill].dG1       = (double)dGx * 0.01;
        storage[storage_fill].dG2       = (double)dGy * 0.01;
        storage[storage_fill].energy    = (double)inter * 0.01;
        storage[storage_fill].structure = struc;
        storage[storage_fill].inactive  = 0;
        storage[storage_fill].processed = 0;

        storage_fill++;

        if (storage_fill == storage_size - 1) {
          storage_size *= 1.4;
          storage       = (vrna_pkplex_t *)vrna_realloc(storage,
                                                        sizeof(vrna_pkplex_t) *
                                                        storage_size);
        }
      }
    }
  }

  free_array(c3, n, max_interaction_length);

  /* resize to space actually required */
  if (storage_fill > 0) {
    storage = (vrna_pkplex_t *)vrna_realloc(storage,
                                            sizeof(vrna_pkplex_t) *
                                            (storage_fill + 1));

    /* add end-of-list identifier */
    storage[storage_fill].structure = NULL;
    storage[storage_fill].inactive  = 1;
  } else {
    free(storage);
    storage = NULL;
  }

  return storage;
}


PRIVATE char *
backtrack_XS(vrna_fold_compound_t *fc,
             int        k,
             int        l,
             const int  i,
             const int  j,
             const int  max_interaction_length,
             int        ***c3)
{
  /* backtrack structure going backwards from i, and forwards from j
   * return structure in bracket notation with & as separator */
  short         *S, *S1;
  int           p, q, type, type2, E, traced, i0, j0, *rtype;
  char          *st1, *st2, *struc;
  vrna_param_t  *P;
  vrna_md_t     *md;

  S               = fc->sequence_encoding2;
  S1              = fc->sequence_encoding;
  P               = fc->params;
  md              = &(P->model_details);
  rtype           = &(md->rtype[0]);
  st1             = (char *)vrna_alloc(sizeof(char) * (i - k + 2));
  st1[i - k + 1]  = '\0';
  st2             = (char *)vrna_alloc(sizeof(char) * (l - j + 2));
  st2[l - j + 1]  = '\0';

  i0  = k;
  j0  = l;
  while (k <= i && l >= j) {
    E           = c3[j - 1][max_interaction_length - i + k - 1][l - j];
    traced      = 0;
    st1[k - i0] = '(';
    st2[l - j]  = ')';

    type = md->pair[S[k]][S[l]];

    if (!type)
      vrna_message_error("backtrack failed in fold duplex bli");

    for (p = k + 1; p <= i; p++) {
      for (q = l - 1; q >= j; q--) {
        int LE;
        if (p - k + l - q - 2 > MAXLOOP)
          break;

        type2 = md->pair[S[q]][S[p]];

        if (!type2)
          continue;

        LE = E_IntLoop(p - k - 1,
                       l - q - 1,
                       type,
                       type2,
                       S1[k + 1],
                       S1[l - 1],
                       S1[p - 1],
                       S1[q + 1],
                       P);
        if (E == c3[j - 1][max_interaction_length - i + p - 1][q - j] + LE) {
          traced  = 1;
          k       = p;
          l       = q;
          break;
        }
      }
      if (traced)
        break;
    }
    if (!traced) {
      E -= vrna_E_ext_stem(rtype[type],
                           S1[l - 1],
                           S1[k + 1],
                           P);

      if (E != 0)
        vrna_message_error("backtrack failed in fold duplex bal");
      else
        break;
    }
  }
  struc = (char *)vrna_alloc(k - i0 + 1 + j0 - l + 1 + 2);

  for (p = 0; p <= i - i0; p++)
    if (!st1[p])
      st1[p] = '.';

  for (p = 0; p <= j0 - j; p++)
    if (!st2[p])
      st2[p] = '.';

  strcpy(struc, st1);
  strcat(struc, "&");
  strcat(struc, st2);
  free(st1);
  free(st2);
  return struc;
}
