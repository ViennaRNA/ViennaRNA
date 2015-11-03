/*
 * Example implementation for including ligand binding to hairpins, or
 * interior loops, with known sequence motif and binding free energy
 * utilizing generic soft constraint feature
 *
 * (c) 2015 Ronny Lorenz - ViennaRNA Package
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/ligand.h"


typedef struct{
  int i;
  int j;
  int k;
  int l;
} quadruple_position;

typedef struct{
  char  *seq_motif_5;
  char  *seq_motif_3;
  char  *struct_motif_5;
  char  *struct_motif_3;
  int   energy;
  int   energy_alt;

  quadruple_position *positions;
} ligand_data;

static void
delete_ligand_data(void *data){

  ligand_data *ldata = (ligand_data *)data;

  free(ldata->seq_motif_5);
  free(ldata->seq_motif_3);
  free(ldata->struct_motif_5);
  free(ldata->struct_motif_3);
  free(ldata->positions);

  free(data);
}

static void
scanForIntMotif(char *seq, char *motif1, char *motif2, quadruple_position **pos){

  int   i, j, k, l, l1, l2, n, cnt, cnt2;
  char  *ptr;
  
  n     = (int) strlen(seq);
  l1    = (int) strlen(motif1);
  l2    = (int) strlen(motif2);
  cnt   = 0;
  cnt2  = 5; /* initial guess how many matching motifs we might encounter */

  *pos = (quadruple_position *)vrna_alloc(sizeof(quadruple_position) * cnt2);

  for(i = 0; i <= n - l1 - l2; i++){
    if(seq[i] == motif1[0]){
      for(j = i+1; j < i + l1; j++){
        if(seq[j] == motif1[j-i]){
          continue;
        }
        else goto next_i;
      }
      for(k = j + 1; k <= n - l2; k++){
        if(seq[k] == motif2[0]){
          for(l = k + 1; l < k + l2; l++){
            if(seq[l] == motif2[l-k]){
              continue;
            }
            else goto next_k;
          }
          /* we found a quadruple, so store it */
          (*pos)[cnt].i   = i + 1;
          (*pos)[cnt].j   = l;
          (*pos)[cnt].k   = j;
          (*pos)[cnt++].l = k + 1;

          /* allocate more memory if necessary */
          if(cnt == cnt2){
            cnt2 *= 2;
            *pos = (quadruple_position *)vrna_realloc(*pos, sizeof(quadruple_position) * cnt2);
          }
        }
/* early exit from l loop */
next_k: continue;
      }
    }
/* early exit from j loop */
next_i: continue;
  }

  /* reallocate to actual size */
  *pos = (quadruple_position *)vrna_realloc(*pos, sizeof(quadruple_position) * (cnt + 1));

  /* set end marker */
  (*pos)[cnt].i = (*pos)[cnt].j = (*pos)[cnt].k = (*pos)[cnt++].l = 0;
}

static int
AptamerContrib(int i, int j, int k, int l, char d, void *data){

  quadruple_position  *pos;
  ligand_data         *ldata;

  if(d & VRNA_DECOMP_PAIR_IL){
    ldata = (ligand_data *)data;
    for(pos = ((ligand_data *)data)->positions; pos->i; pos++){
      if((pos->i == i) && (pos->j == j) && (pos->k == k) && (pos->l == l)){
        return ldata->energy;
      }
    }
  }

  return 0;
}

static FLT_OR_DBL
expAptamerContrib(int i, int j, int k, int l, char d, void *data){

  quadruple_position  *pos;
  ligand_data         *ldata;
  FLT_OR_DBL          exp_e;
  double              kT;


  if(d & VRNA_DECOMP_PAIR_IL){
    ldata = (ligand_data *)data;
    exp_e = 1.;
    kT    = (37. + K0) * GASCONST;

    for(pos = ldata->positions; pos->i; pos++){
      if((pos->i == i) && (pos->j == j) && (pos->k == k) && (pos->l == l)){
        exp_e =   exp(-ldata->energy*10./kT);
        exp_e +=  exp(-ldata->energy_alt*10./kT); /* add alternative, i.e. unbound ligand */
      }
    }
  }

  return exp_e;
}

PUBLIC int
vrna_sc_add_hi_motif( vrna_fold_compound_t *vc,
                      const char *seq,
                      const char *structure,
                      double energy,
                      unsigned int options){

    int                   i, cp, cp2;
    char                  *sequence, *motif, *motif2;
    vrna_md_t             md;
    ligand_data           *ldata;
    vrna_fold_compound_t  *tmp_vc;

    set_model_details(&md);

    ldata                 = vrna_alloc(sizeof(ligand_data));
    ldata->struct_motif_5 = NULL;
    ldata->struct_motif_3 = NULL;
    ldata->positions      = NULL;

    sequence  = vrna_cut_point_remove(seq, &cp);        /* ligand sequence motif  */
    motif     = vrna_cut_point_remove(structure, &cp2); /* ligand structure motif */

    if(cp != cp2){
      vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: Cutpoint in sequence and structure motif differ!");
      return 0;
    } else if(strlen(seq) != strlen(structure)){
      vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: length of sequence and structure motif differ!");
      return 0;
    }

    if(cp > 0){
      ldata->seq_motif_5 = vrna_alloc(sizeof(char) * cp);
      strncpy(ldata->seq_motif_5, sequence, cp - 1);
      ldata->seq_motif_5[cp - 1] = '\0';
      ldata->struct_motif_5 = vrna_alloc(sizeof(char) * cp);
      strncpy(ldata->struct_motif_5, motif, cp - 1);
      ldata->struct_motif_5[cp - 1] = '\0';

      ldata->seq_motif_3 = vrna_alloc(sizeof(char) * (strlen(sequence) - cp + 2));
      strncpy(ldata->seq_motif_3, sequence + cp - 1, (strlen(sequence) - cp + 1));
      ldata->seq_motif_3[strlen(sequence) - cp + 1] = '\0';
      ldata->struct_motif_3 = vrna_alloc(sizeof(char) * (strlen(motif) - cp + 2));
      strncpy(ldata->struct_motif_3, motif + cp - 1, (strlen(motif) - cp + 1));
      ldata->struct_motif_3[strlen(motif) - cp + 1] = '\0';

      /* compute alternative contribution and energy correction */
      motif2  = vrna_alloc(sizeof(char) * strlen(structure));
      memset(motif2, '.', strlen(structure) - 1);
      motif2[strlen(structure)] = '\0';
      /* insert corresponding interior loop motif (....(&)...) */
      motif2[0] = '(';
      motif2[cp-2] = '(';
      motif2[cp-1] = ')';
      motif2[strlen(structure) - 2] = ')';
      motif2[strlen(structure) - 1] = '\0';

      /* create temporary vrna_fold_compound for energy evaluation */
      tmp_vc = vrna_fold_compound(seq, &md, VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY);
      int alt     = (int)(vrna_eval_structure(tmp_vc, motif2) * 100.);
      int corr    = (int)(vrna_eval_structure(tmp_vc, motif) * 100.);
      energy     += (double)(corr - alt) / 100.;

      ldata->energy     = (int)(energy*100.);
      ldata->energy_alt = alt;
    } else {
      ldata->seq_motif_5 = sequence;
      ldata->seq_motif_3 = NULL;
    }

    /* scan for motif positions */
    scanForIntMotif(vc->sequence, ldata->seq_motif_5, ldata->seq_motif_3, &(ldata->positions));

    /* add generalized soft-constraints to predict aptamer-ligand complexes */
    vrna_sc_add_data(vc, (void *)ldata, &delete_ligand_data);
    if(options & VRNA_OPTION_MFE)
      vrna_sc_add_f(vc, &AptamerContrib);
    if(options & VRNA_OPTION_PF)
      vrna_sc_add_exp_f(vc, &expAptamerContrib);

    return 1; /* success */
}
