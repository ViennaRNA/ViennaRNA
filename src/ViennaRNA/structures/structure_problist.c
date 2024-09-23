/*
 *  ViennaRNA/utils/structures.c
 *
 *  Various functions to convert, parse, encode secondary structures
 *
 *  c  Ivo L Hofacker, Walter Fontana, Ronny Lorenz
 *              Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/partfunc/gquad.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/structures/problist.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE vrna_ep_t *
wrap_get_plist(vrna_mx_pf_t     *matrices,
               int              length,
               int              *index,
               short            *S,
               vrna_exp_param_t *pf_params,
               double           cut_off);


PRIVATE vrna_ep_t *
wrap_plist(vrna_fold_compound_t *vc,
           double               cut_off);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC vrna_ep_t *
vrna_plist(const char *struc,
           float      pr)
{
  /* convert bracket string to plist */
  short     *pt;
  int       i, k = 0, size, n;
  vrna_ep_t *gpl, *ptr, *pl;

  pl = NULL;

  if (struc) {
    size  = strlen(struc);
    n     = 2;

    pt  = vrna_ptable(struc);
    pl  = (vrna_ep_t *)vrna_alloc(n * size * sizeof(vrna_ep_t));
    for (i = 1; i < size; i++) {
      if (pt[i] > i) {
        (pl)[k].i       = i;
        (pl)[k].j       = pt[i];
        (pl)[k].p       = pr;
        (pl)[k++].type  = VRNA_PLIST_TYPE_BASEPAIR;
      }
    }

    gpl = get_plist_gquad_from_db(struc, pr);
    for (ptr = gpl; ptr->i != 0; ptr++) {
      if (k == n * size - 1) {
        n   *= 2;
        pl  = (vrna_ep_t *)vrna_realloc(pl, n * size * sizeof(vrna_ep_t));
      }

      (pl)[k].i       = ptr->i;
      (pl)[k].j       = ptr->j;
      (pl)[k].p       = ptr->p;
      (pl)[k++].type  = ptr->type;
    }
    free(gpl);

    (pl)[k].i       = 0;
    (pl)[k].j       = 0;
    (pl)[k].p       = 0.;
    (pl)[k++].type  = 0.;
    free(pt);
    pl = (vrna_ep_t *)vrna_realloc(pl, k * sizeof(vrna_ep_t));
  }

  return pl;
}


PUBLIC vrna_ep_t *
vrna_plist_from_probs(vrna_fold_compound_t  *vc,
                      double                cut_off)
{
  if (!vc)
    vrna_log_warning("vrna_pl_get_from_pr: run vrna_pf_fold first!");
  else if (!vc->exp_matrices->probs)
    vrna_log_warning("vrna_pl_get_from_pr: probs==NULL!");
  else
    return wrap_plist(vc, cut_off);

  return NULL;
}


PUBLIC int
vrna_plist_append(vrna_ep_t       **target,
                  const vrna_ep_t *list)
{
  int             size1, size2;
  const vrna_ep_t *ptr;

  if ((target) && (list)) {
    size1 = size2 = 0;

    if (*target)
      for (ptr = *target; ptr->i; size1++, ptr++);

    for (ptr = list; ptr->i; size2++, ptr++);

    *target = (vrna_ep_t *)vrna_realloc(*target, sizeof(vrna_ep_t) * (size1 + size2 + 1));

    if (*target) {
      memcpy(*target + size1, list, sizeof(vrna_ep_t) * size2);
      (*target)[size1 + size2].i    = (*target)[size1 + size2].j = 0;
      (*target)[size1 + size2].type = 0;
      return 1;
    }
  }

  return 0;
}


PUBLIC plist *
get_plist_gquad_from_db(const char  *structure,
                        float       pr)
{
  unsigned int  x, L, l[3], n, size, ge, ee, actual_size, gb;
  plist         *pl;

  actual_size = 0;
  ge          = 0;
  n           = 2;
  size        = strlen(structure);
  pl          = (plist *)vrna_alloc(n * size * sizeof(plist));

  while ((ee = vrna_gq_parse(structure + ge, &L, l)) > 0) {
    ge  += ee;

    if (4 * L + l[0] + l[1] + l[2] > ee) {
      gb = size + ge - L * 4 - l[0] - l[1] - l[2] + 1;
    } else {
      gb = ge - L * 4 - l[0] - l[1] - l[2] + 1;
    }

    /* add pseudo-base pair enclosing gquad */
    if (actual_size >= n * size - 5) {
      n   *= 2;
      pl  = (plist *)vrna_realloc(pl, n * size * sizeof(plist));
    }

    pl[actual_size].i       = gb;
    pl[actual_size].j       = ge;
    pl[actual_size].p       = pr;
    pl[actual_size++].type  = VRNA_PLIST_TYPE_GQUAD;

    for (x = 0; x < L; x++) {
      if (actual_size >= n * size - 5) {
        n   *= 2;
        pl  = (plist *)vrna_realloc(pl, n * size * sizeof(plist));
      }

      pl[actual_size].i       = (gb + x - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - L + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (gb + x - 1) % (size) + 1;
      pl[actual_size].j       = (gb + x + l[0] + L - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (gb + x + l[0] + L - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - 2 * L - l[2] + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;

      pl[actual_size].i       = (ge + x - 2 * L - l[2] + 1 - 1) % (size) + 1;
      pl[actual_size].j       = (ge + x - L + 1 - 1) % (size) + 1;
      pl[actual_size].p       = pr;
      pl[actual_size++].type  = VRNA_PLIST_TYPE_TRIPLE;
    }
  }

  pl[actual_size].i   = pl[actual_size].j = 0;
  pl[actual_size++].p = 0;
  pl                  = (plist *)vrna_realloc(pl, actual_size * sizeof(plist));
  return pl;
}


/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE vrna_ep_t *
wrap_get_plist(vrna_mx_pf_t     *matrices,
               int              length,
               int              *index,
               short            *S,
               vrna_exp_param_t *pf_params,
               double           cut_off)
{
  int         i, j, k, n, count, gquad;
  FLT_OR_DBL  *probs, *scale;
  vrna_ep_t   *pl;
  vrna_smx_csr(FLT_OR_DBL) *q_gq;

  probs = matrices->probs;
  q_gq  = matrices->q_gq;
  scale = matrices->scale;
  gquad = pf_params->model_details.gquad;

  count = 0;
  n     = 2;

  /* first guess of the size needed for pl */
  pl = (vrna_ep_t *)vrna_alloc(n * length * sizeof(vrna_ep_t));

  for (i = 1; i < length; i++) {
    for (j = i + 1; j <= length; j++) {
      /* skip all entries below the cutoff */
      if (probs[index[i] - j] < (FLT_OR_DBL)cut_off)
        continue;

      /* do we need to allocate more memory? */
      if (count == n * length - 1) {
        n   *= 2;
        pl  = (vrna_ep_t *)vrna_realloc(pl, n * length * sizeof(vrna_ep_t));
      }

      /* check for presence of gquadruplex */
      if (gquad && (S[i] == 3) && (S[j] == 3)) {
        /* add probability of a gquadruplex at position (i,j)
         * for dot_plot
         */
        (pl)[count].i       = i;
        (pl)[count].j       = j;
        (pl)[count].p       = (float)probs[index[i] - j];
        (pl)[count++].type  = VRNA_PLIST_TYPE_GQUAD;
        /* now add the probabilies of it's actual pairing patterns */
        vrna_ep_t *inner, *ptr;
        inner = get_plist_gquad_from_pr(S, i, j, q_gq, probs, scale, pf_params);
        for (ptr = inner; ptr->i != 0; ptr++) {
          if (count == n * length - 1) {
            n   *= 2;
            pl  = (vrna_ep_t *)vrna_realloc(pl, n * length * sizeof(vrna_ep_t));
          }

          /* check if we've already seen this pair */
          for (k = 0; k < count; k++)
            if (((pl)[k].i == ptr->i) &&
                ((pl)[k].j == ptr->j))
              break;

          (pl)[k].i     = ptr->i;
          (pl)[k].j     = ptr->j;
          (pl)[k].type  = ptr->type;
          if (k == count) {
            (pl)[k].p = ptr->p;
            count++;
          } else {
            (pl)[k].p += ptr->p;
          }
        }
      } else {
        (pl)[count].i       = i;
        (pl)[count].j       = j;
        (pl)[count].p       = (float)probs[index[i] - j];
        (pl)[count++].type  = VRNA_PLIST_TYPE_BASEPAIR;
      }
    }
  }
  /* mark the end of pl */
  (pl)[count].i     = 0;
  (pl)[count].j     = 0;
  (pl)[count].type  = 0;
  (pl)[count++].p   = 0.;
  /* shrink memory to actual size needed */
  pl = (vrna_ep_t *)vrna_realloc(pl, count * sizeof(vrna_ep_t));

  return pl;
}


PRIVATE vrna_ep_t *
wrap_plist(vrna_fold_compound_t *vc,
           double               cut_off)
{
  short             *S;
  unsigned int      with_gquad, with_circ, i, j, k, n, m, length, count;
  int               *index;
  FLT_OR_DBL        *probs;
  vrna_ep_t         *pl;
  vrna_mx_pf_t      *matrices;
  vrna_exp_param_t  *pf_params;
  vrna_smx_csr(FLT_OR_DBL)  *p_gq;

  S           = (vc->type == VRNA_FC_TYPE_SINGLE) ? vc->sequence_encoding2 : vc->S_cons;
  index       = vc->iindx;
  length      = vc->length;
  pf_params   = vc->exp_params;
  matrices    = vc->exp_matrices;
  probs       = matrices->probs;
  with_gquad  = pf_params->model_details.gquad;
  with_circ   = pf_params->model_details.circ;

  count = 0;
  n     = 2;

  /* first guess of the size needed for pl */
  pl = (vrna_ep_t *)vrna_alloc(n * length * sizeof(vrna_ep_t));

  for (i = 1; i < length; i++) {
    for (j = i + 1; j <= length; j++) {
      /* skip all entries below the cutoff */
      if (probs[index[i] - j] < (FLT_OR_DBL)cut_off)
        continue;

      /* do we need to allocate more memory? */
      if (count == n * length - 1) {
        n   *= 2;
        pl  = (vrna_ep_t *)vrna_realloc(pl, n * length * sizeof(vrna_ep_t));
      }

      /* check for presence of gquadruplex */
      if ((with_gquad) &&
          (S[i] == 3) &&
          (S[j] == 3)) {
        /* add probability of a gquadruplex at position (i,j)
         * for dot_plot
         */
        (pl)[count].i       = i;
        (pl)[count].j       = j;
        (pl)[count].p       = (float)probs[index[i] - j];
        (pl)[count++].type  = VRNA_PLIST_TYPE_GQUAD;

        /* now add the probabilies of it's actual pairing patterns */
        vrna_ep_t *inner, *ptr;
        inner = vrna_plist_gquad_from_pr(vc, i, j);
        for (ptr = inner; ptr->i != 0; ptr++) {
          if (count == n * length - 1) {
            n   *= 2;
            pl  = (vrna_ep_t *)vrna_realloc(pl, n * length * sizeof(vrna_ep_t));
          }

          /* check if we've already seen this pair */
          for (k = 0; k < count; k++)
            if (((pl)[k].i == ptr->i) &&
                ((pl)[k].j == ptr->j) &&
                ((pl)[k].type == VRNA_PLIST_TYPE_TRIPLE))
              break;

          (pl)[k].i     = ptr->i;
          (pl)[k].j     = ptr->j;
          (pl)[k].type  = ptr->type;
          if (k == count) {
            (pl)[k].p = ptr->p;
            count++;
          } else {
            (pl)[k].p += ptr->p;
          }
        }
        free(inner);
      } else {
        (pl)[count].i       = i;
        (pl)[count].j       = j;
        (pl)[count].p       = (float)probs[index[i] - j];
        (pl)[count++].type  = VRNA_PLIST_TYPE_BASEPAIR;
      }
    }
  }

  if ((with_gquad) &&
      (with_circ) &&
      (vc->exp_matrices->p_gq)) {
    p_gq = vc->exp_matrices->p_gq;
    unsigned int imin = 2;
    if (imin + VRNA_GQUAD_MAX_BOX_SIZE - 1 <= length)
       imin = length - VRNA_GQUAD_MAX_BOX_SIZE + 1;

    for (i = imin; i <= length; i++) {
      unsigned int jmin = 1;
      unsigned int jmax = VRNA_GQUAD_MAX_BOX_SIZE - 1 - (length - i);
      if (jmax >= length)
        jmax = i - 1;

      for (j = jmin; j <= jmax; j++) {
        FLT_OR_DBL p_g;
#ifndef VRNA_DISABLE_C11_FEATURES
        if ((p_g = vrna_smx_csr_get(p_gq, i, j, 0.)) >= (FLT_OR_DBL)cut_off) {
#else
        if ((p_g = vrna_smx_csr_FLT_OR_DBL_get(p_gq, i, j, 0.)) >= (FLT_OR_DBL)cut_off) {
#endif
          (pl)[count].i       = i;
          (pl)[count].j       = j;
          (pl)[count].p       = (float)p_g;
          (pl)[count++].type  = VRNA_PLIST_TYPE_GQUAD;

          /* now add the probabilies of it's actual pairing patterns */
          vrna_ep_t *inner, *ptr;
          inner = vrna_plist_gquad_from_pr(vc, i, j);
          for (ptr = inner; ptr->i != 0; ptr++) {
            if (count == n * length - 1) {
              n   *= 2;
              pl  = (vrna_ep_t *)vrna_realloc(pl, n * length * sizeof(vrna_ep_t));
            }

            /* check if we've already seen this pair */
            for (k = 0; k < count; k++)
              if (((pl)[k].i == ptr->i) &&
                  ((pl)[k].j == ptr->j) &&
                  ((pl)[k].type == VRNA_PLIST_TYPE_TRIPLE))
                break;

            (pl)[k].i     = ptr->i;
            (pl)[k].j     = ptr->j;
            (pl)[k].type  = ptr->type;
            if (k == count) {
              (pl)[k].p = ptr->p;
              count++;
            } else {
              (pl)[k].p += ptr->p;
            }
          }
          free(inner);
        }
      }
    }
  }

  /* check unstructured domains */
  if (vc->domains_up) {
    vrna_ud_t *domains_up;
    domains_up = vc->domains_up;

    if (domains_up->probs_get) {
      for (i = 1; i <= length; i++)
        for (m = 0; m < domains_up->motif_count; m++) {
          FLT_OR_DBL pp;
          j   = i + domains_up->motif_size[m] - 1;
          pp  = 0.;
          pp  += domains_up->probs_get(vc,
                                       i,
                                       j,
                                       VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                                       m,
                                       domains_up->data);
          pp += domains_up->probs_get(vc,
                                      i,
                                      j,
                                      VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                      m,
                                      domains_up->data);
          pp += domains_up->probs_get(vc,
                                      i,
                                      j,
                                      VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                      m,
                                      domains_up->data);
          pp += domains_up->probs_get(vc,
                                      i,
                                      j,
                                      VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                      m,
                                      domains_up->data);
          if (pp >= (FLT_OR_DBL)cut_off) {
            /* do we need to allocate more memory? */
            if (count == n * length - 1) {
              n   *= 2;
              pl  = (vrna_ep_t *)vrna_realloc(pl, n * length * sizeof(vrna_ep_t));
            }

            (pl)[count].i       = i;
            (pl)[count].j       = j;
            (pl)[count].p       = (float)pp;
            (pl)[count++].type  = VRNA_PLIST_TYPE_UD_MOTIF;
          }
        }
    }
  }

  /* mark the end of pl */
  (pl)[count].i     = 0;
  (pl)[count].j     = 0;
  (pl)[count].type  = 0;
  (pl)[count++].p   = 0.;
  /* shrink memory to actual size needed */
  pl = (vrna_ep_t *)vrna_realloc(pl, count * sizeof(vrna_ep_t));

  return pl;
}


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

PUBLIC void
assign_plist_from_pr(vrna_ep_t  **pl,
                     FLT_OR_DBL *probs,
                     int        length,
                     double     cut_off)
{
  int               *index;
  vrna_mx_pf_t      *matrices;
  vrna_md_t         md;
  vrna_exp_param_t  *pf_params;

  index     = vrna_idx_row_wise(length);
  matrices  = (vrna_mx_pf_t *)vrna_alloc(sizeof(vrna_mx_pf_t));

  set_model_details(&md);
  md.gquad        = 0;
  pf_params       = vrna_exp_params(&md);
  matrices->probs = probs;

  *pl = wrap_get_plist(matrices,
                       length,
                       index,
                       NULL,
                       pf_params,
                       cut_off);

  free(index);
  free(pf_params);
  free(matrices);
}


PUBLIC void
assign_plist_from_db(vrna_ep_t  **pl,
                     const char *struc,
                     float      pr)
{
  *pl = vrna_plist(struc, pr);
}

#endif
