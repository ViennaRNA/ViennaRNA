/*
 *                centroid structure prediction
 *
 *                Ivo L Hofacker + Ronny Lorenz
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/partfunc/gquad.h"
#include "ViennaRNA/structures/centroid.h"

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

/* compute the centroid structure of the ensemble, i.e. the strutcure
 * with the minimal average distance to all other structures
 * <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
 * Thus, the centroid is simply the structure containing all pairs with
 * p_ij>0.5
 */
PUBLIC char *
vrna_centroid_from_plist(int        length,
                         double     *dist,
                         vrna_ep_t  *pl)
{
  int   i;
  char  *centroid;

  if (pl == NULL) {
    vrna_log_warning("vrna_centroid_from_plist: "
                     "pl == NULL!");
    return NULL;
  }

  *dist     = 0.;
  centroid  = (char *)vrna_alloc((length + 1) * sizeof(char));
  for (i = 0; i < length; i++)
    centroid[i] = '.';
  for (i = 0; pl[i].i > 0; i++) {
    if ((pl[i].p) > 0.5) {
      centroid[pl[i].i - 1] = '(';
      centroid[pl[i].j - 1] = ')';
      *dist                 += (1 - pl[i].p);
    } else {
      *dist += pl[i].p;
    }
  }
  centroid[length] = '\0';
  return centroid;
}


/* compute the centroid structure of the ensemble, i.e. the strutcure
 * with the minimal average distance to all other structures
 * <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
 * Thus, the centroid is simply the structure containing all pairs with
 * p_ij>0.5
 */
PUBLIC char *
vrna_centroid_from_probs(int        length,
                         double     *dist,
                         FLT_OR_DBL *probs)
{
  int         i, j;
  FLT_OR_DBL  p;
  char        *centroid;
  int         *index = vrna_idx_row_wise(length);

  if (probs == NULL) {
    vrna_log_warning("vrna_centroid_from_probs: "
                     "probs == NULL!");
    return NULL;
  }

  *dist     = 0.;
  centroid  = (char *)vrna_alloc((length + 1) * sizeof(char));
  for (i = 0; i < length; i++)
    centroid[i] = '.';
  for (i = 1; i <= length; i++)
    for (j = i + 1; j <= length; j++) {
      if ((p = probs[index[i] - j]) > 0.5) {
        centroid[i - 1] = '(';
        centroid[j - 1] = ')';
        *dist           += (1 - p);
      } else {
        *dist += p;
      }
    }
  free(index);
  centroid[length] = '\0';
  return centroid;
}


/* compute the centroid structure of the ensemble, i.e. the strutcure
 * with the minimal average distance to all other structures
 * <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
 * Thus, the centroid is simply the structure containing all pairs with
 * p_ij>0.5
 */
PUBLIC char *
vrna_centroid(vrna_fold_compound_t  *fc,
              double                *dist)
{
  char              *centroid;
  short             *S;
  unsigned int      L, l[3], i, j, k, kmax, n;
  int               *my_iindx;
  FLT_OR_DBL        p, *probs;
  vrna_mx_pf_t      *matrices;
  vrna_exp_param_t  *pf_params;
  vrna_md_t         *md;

  if (!fc) {
    vrna_log_warning("vrna_fold_compound_t missing!");
    return NULL;
  } else if (!dist) {
    vrna_log_error("pointer to centroid distance is missing");
    return NULL;
  } else if (!fc->exp_matrices->probs) {
    vrna_log_warning("probs == NULL!");
    return NULL;
  }

  n         = fc->length;
  pf_params = fc->exp_params;
  md        = &(pf_params->model_details);
  S         = (fc->type == VRNA_FC_TYPE_SINGLE) ? fc->sequence_encoding2 : fc->S_cons;
  my_iindx  = fc->iindx;
  matrices  = fc->exp_matrices;
  probs     = matrices->probs;
  *dist     = 0.;

  centroid = (char *)vrna_alloc((n + 1) * sizeof(char));

  for (i = 0; i < n; i++)
    centroid[i] = '.';
  for (i = 1; i <= n; i++)
    for (j = i + 1; j <= n; j++) {
      if ((p = probs[my_iindx[i] - j]) > 0.5) {
        if (md->gquad) {
          /* check for presence of gquadruplex */
          if ((S[i] == 3) && (S[j] == 3)) {
            vrna_get_gquad_pattern_pf(fc, i, j, &L, l);

            if (L > 0)
              vrna_db_insert_gq(centroid, i, L, l, n);
            else
              vrna_log_error("failed to detect G-Quadruplex pattern");

            /* skip everything within the gquad */
            i     = j;
            j     = j + 1;
            *dist += (1 - p); /* right? */
            break;
          }
        }

        /* regular base pair */
        centroid[i - 1] = '(';
        centroid[j - 1] = ')';
        *dist           += (1 - p);
      } else {
        *dist += p;
      }
    }

  if ((md->circ) &&
      (md->gquad) &&
      (matrices->p_gq)) {
#ifndef VRNA_DISABLE_C11_FEATURES
    kmax = vrna_smx_csr_get_size(matrices->p_gq);
#else
    kmax = vrna_smx_csr_FLT_OR_DBL_get_size(matrices->p_gq);
#endif

    /* add gquads that wrap around the n,1 junction */
    for (k = 0; k < kmax; k++) {
#ifndef VRNA_DISABLE_C11_FEATURES
      p = vrna_smx_csr_get_entry(matrices->p_gq, k, &i, &j, 0.);
#else
      p = vrna_smx_csr_FLT_OR_DBL_get_entry(matrices->p_gq, k, &i, &j, 0.);
#endif
      if (p > 0.5) {
        /* get the gquad pattern */
        vrna_get_gquad_pattern_pf(fc, i, j, &L, l);

        if (L > 0)
          vrna_db_insert_gq(centroid, i, L, l, n);
        else
          vrna_log_error("failed to detect G-Quadruplex pattern");

        *dist += (1 - p);
      }
    }
  }

  centroid[n] = '\0';
  return centroid;
}


/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */


/* this function is a threadsafe replacement for centroid() */
PUBLIC char *
get_centroid_struct_pl(int        length,
                       double     *dist,
                       vrna_ep_t  *pl)
{
  return vrna_centroid_from_plist(length, dist, pl);
}


/* this function is a threadsafe replacement for centroid() */
PUBLIC char *
get_centroid_struct_pr(int        length,
                       double     *dist,
                       FLT_OR_DBL *probs)
{
  return vrna_centroid_from_probs(length, dist, probs);
}
