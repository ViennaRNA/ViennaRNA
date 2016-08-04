/*
                  centroid structure prediction

                  Ivo L Hofacker + Ronny Lorenz
                  Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/centroid.h"

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

PUBLIC char *
vrna_centroid_from_plist( int length,
                          double *dist,
                          vrna_plist_t *pl){

  /* compute the centroid structure of the ensemble, i.e. the strutcure
     with the minimal average distance to all other structures
     <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
     Thus, the centroid is simply the structure containing all pairs with
     p_ij>0.5 */
  int i;
  char *centroid;

  if (pl==NULL)
    vrna_message_error("vrna_centroid_from_plist@centroid.c: pl==NULL!");

  *dist = 0.;
  centroid = (char *) vrna_alloc((length+1)*sizeof(char));
  for (i=0; i<length; i++) centroid[i]='.';
  for (i=0; pl[i].i>0; i++){
    if ((pl[i].p)>0.5) {
      centroid[pl[i].i-1] = '(';
      centroid[pl[i].j-1] = ')';
      *dist += (1-pl[i].p);
    } else
      *dist += pl[i].p;
  }
  centroid[length] = '\0';
  return centroid;
}

PUBLIC char *
vrna_centroid_from_probs( int length,
                          double *dist,
                          FLT_OR_DBL *probs){

  /* compute the centroid structure of the ensemble, i.e. the strutcure
     with the minimal average distance to all other structures
     <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
     Thus, the centroid is simply the structure containing all pairs with
     p_ij>0.5 */
  int i,j;
  FLT_OR_DBL p;
  char  *centroid;
  int   *index = vrna_idx_row_wise(length);

  if (probs == NULL)
    vrna_message_error("vrna_centroid_from_probs@centroid.c: probs==NULL!");

  *dist = 0.;
  centroid = (char *) vrna_alloc((length+1)*sizeof(char));
  for (i=0; i<length; i++) centroid[i]='.';
  for (i=1; i<=length; i++)
    for (j=i+TURN+1; j<=length; j++) {
      if ((p=probs[index[i]-j])>0.5) {
        centroid[i-1] = '(';
        centroid[j-1] = ')';
        *dist += (1-p);
      } else
        *dist += p;
    }
  free(index);
  centroid[length] = '\0';
  return centroid;
}

PUBLIC char *
vrna_centroid(vrna_fold_compound_t *vc,
             double *dist){

  /* compute the centroid structure of the ensemble, i.e. the strutcure
     with the minimal average distance to all other structures
     <d(S)> = \sum_{(i,j) \in S} (1-p_{ij}) + \sum_{(i,j) \notin S} p_{ij}
     Thus, the centroid is simply the structure containing all pairs with
     p_ij>0.5 */
  int i,j, k, length;
  FLT_OR_DBL p;
  char  *centroid;
  short *S;
  vrna_mx_pf_t      *matrices;
  FLT_OR_DBL        *probs;
  int               *my_iindx;
  vrna_exp_param_t  *pf_params;


  if(!vc){
    vrna_message_error("vrna_centroid@centroid.c: run vrna_pf_fold first!");
  } else if( !vc->exp_matrices->probs){
    vrna_message_error("vrna_centroid@centroid.c: probs==NULL!");
  }

  length      = vc->length;
  pf_params   = vc->exp_params;
  S           = vc->sequence_encoding2;
  my_iindx    = vc->iindx;

  matrices    = vc->exp_matrices;
  probs       = matrices->probs;

  *dist = 0.;
  centroid = (char *) vrna_alloc((length+1)*sizeof(char));
  for (i=0; i<length; i++) centroid[i]='.';
  for (i=1; i<=length; i++)
    for (j=i+TURN+1; j<=length; j++) {
      if ((p=probs[my_iindx[i]-j])>0.5) {
        if(pf_params->model_details.gquad){
          /* check for presence of gquadruplex */
          if((S[i] == 3) && (S[j] == 3)){
            int L, l[3];
            get_gquad_pattern_pf(S, i, j, pf_params, &L, l);
            for(k=0;k<L;k++){
              centroid[i+k-1]\
              = centroid[i+k+L+l[0]-1]\
              = centroid[i+k+2*L+l[0]+l[1]-1]\
              = centroid[i+k+3*L+l[0]+l[1]+l[2]-1]\
              = '+';
            }
            /* skip everything within the gquad */
            i = j; j = j+TURN+1;
            *dist += (1-p); /* right? */
            break;
          }
        }
        /* regular base pair */
        centroid[i-1] = '(';
        centroid[j-1] = ')';
        *dist += (1-p);
      } else
        *dist += p;
    }

  centroid[length] = '\0';
  return centroid;
}


/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/


/* this function is a threadsafe replacement for centroid() */
PUBLIC char *
get_centroid_struct_pl( int length,
                        double *dist,
                        vrna_plist_t *pl){

  return vrna_centroid_from_plist(length, dist, pl);
}

/* this function is a threadsafe replacement for centroid() */
PUBLIC char *
get_centroid_struct_pr( int length,
                        double *dist,
                        FLT_OR_DBL *probs){

  return vrna_centroid_from_probs(length, dist, probs);
}

