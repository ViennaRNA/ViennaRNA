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
  int   pair_count;
  vrna_basepair_t *pairs;

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
  free(ldata->pairs);

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

static vrna_basepair_t *
backtrack_hi_motif(int i, int j, int k, int l, char d, void *data){

  int                 bp_size = 15;
  vrna_basepair_t     *pairs = NULL;
  quadruple_position  *pos;
  ligand_data         *ldata;

  if(d & VRNA_DECOMP_PAIR_IL){
    ldata = (ligand_data *)data;
    for(pos = ldata->positions; pos->i; pos++){
      if((pos->i == i) && (pos->j == j) && (pos->k == k) && (pos->l == l)){
        /* found motif in our list, lets create pairs */
        char  *ptr;
#if 0
        int   actual_size = 0;
        pairs = vrna_alloc(sizeof(vrna_basepair_t) * bp_size);

        for(ptr=ldata->struct_motif_5; *ptr != '\0'; ptr++, i++){
          if(*ptr == '.'){
            pairs[actual_size].i = pairs[actual_size].j = i;
            actual_size++;
            if(actual_size == bp_size){
              bp_size *= 2;
              pairs = vrna_realloc(pairs, sizeof(vrna_basepair_t) * bp_size);
            }
          }
        }
        for(ptr=ldata->struct_motif_3; *ptr != '\0'; ptr++, l++){
          if(*ptr == '.'){
            pairs[actual_size].i = pairs[actual_size].j = l;
            actual_size++;
            if(actual_size == bp_size){
              bp_size *= 2;
              pairs = vrna_realloc(pairs, sizeof(vrna_basepair_t) * bp_size);
            }
          }
        }
        pairs = vrna_realloc(pairs, sizeof(vrna_basepair_t) * (actual_size + 1));
        pairs[actual_size].i = pairs[actual_size].j = -1;
#else
        pairs = vrna_alloc(sizeof(vrna_basepair_t) * (ldata->pair_count + 1));
        vrna_basepair_t *pptr;
        int             count;
        for(count = 0,pptr = ldata->pairs; pptr->i != -1; pptr++, count++){
          pairs[count].i = i + pptr->i - 1;
          pairs[count].j = (pptr->j < 0) ? j + pptr->j + 1 : i + pptr->j - 1;
        }
        pairs[count].i = pairs[count].j = -1;
#endif

        return pairs;
      }
    }
  }

  return pairs;
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
      motif2  = vrna_alloc(sizeof(char) * (strlen(structure) + 1));
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

      /* pre-construct additional base pairs for backtracking */
      int *stack = vrna_alloc(sizeof(int)*(strlen(seq)+1));
      vrna_basepair_t *pairs = vrna_alloc(sizeof(vrna_basepair_t)*(strlen(seq)+1));
      int pair_count = 0;
      int stack_count = 0;
      for(stack_count = 0, i = 2; i < strlen(ldata->struct_motif_5); i++){ /* go through 5' side of motif */
        if(ldata->struct_motif_5[i-1] == '('){
          stack[stack_count++] = i;
        } else if(ldata->struct_motif_5[i-1] == ')'){
          pairs[pair_count].i = stack[--stack_count];
          pairs[pair_count].j = i;
          pair_count++;
          if(stack_count < 0){
            vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: 5' structure motif contains unbalanced brackets");
            return 0;
          }
        }
      }
      for(i = 2; i < strlen(ldata->struct_motif_3); i++){ /* go through 3' side of motif */
        if(ldata->struct_motif_3[i-1] == '('){
          stack[stack_count++] = -(strlen(ldata->struct_motif_3) - (strlen(ldata->struct_motif_3) - i) - 1);
        } else if(ldata->struct_motif_3[i-1] == ')'){
          pairs[pair_count].i = stack[--stack_count];
          pairs[pair_count].j = -(strlen(ldata->struct_motif_3) - (strlen(ldata->struct_motif_3) - i) - 1);
          pair_count++;
          if(stack_count < 0){
            vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: 3' structure motif contains unbalanced brackets");
            return 0;
          }
        }
      }
      pairs = vrna_realloc(pairs, sizeof(vrna_basepair_t) * (pair_count + 1));
      pairs[pair_count].i = pairs[pair_count].j = -1;
      ldata->pairs = pairs;
      ldata->pair_count = pair_count;
      free(stack);
    } else {
      ldata->seq_motif_5 = sequence;
      ldata->seq_motif_3 = NULL;
    }

    /* scan for motif positions */
    scanForIntMotif(vc->sequence, ldata->seq_motif_5, ldata->seq_motif_3, &(ldata->positions));

    /* add generalized soft-constraints to predict aptamer-ligand complexes */
    vrna_sc_add_data(vc, (void *)ldata, &delete_ligand_data);
    vrna_sc_add_bt(vc, &backtrack_hi_motif);
    if(options & VRNA_OPTION_MFE)
      vrna_sc_add_f(vc, &AptamerContrib);
    if(options & VRNA_OPTION_PF)
      vrna_sc_add_exp_f(vc, &expAptamerContrib);

    return 1; /* success */
}
