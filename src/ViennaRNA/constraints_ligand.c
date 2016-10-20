/*
 * Reference implementation for including ligand binding to hairpins, or
 * interior loops, with known sequence and/or structure motif, and
 * binding free energy utilizing generic soft constraint feature
 *
 * (c) 2015 Ronny Lorenz - ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/constraints.h"
#include "ViennaRNA/eval.h"
#include "ViennaRNA/constraints_ligand.h"


/*
#################################
# PRIVATE DATA STRUCTURES       #
#################################
*/

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

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

static void
split_sequence( const char *string,
                      char **seq1,
                      char **seq2,
                      int cp);

static void
correctMotifContribution( const char *seq,
                          const char *struct_motif,
                          const char *struct_motif_alt,
                          int *contribution,
                          int *contribution_alt,
                          vrna_md_t *md);

static void
delete_ligand_data(void *data);

static int
AptamerContrib(int i, int j, int k, int l, char d, void *data);

static int
AptamerContribHairpin(int i, int j, int k, int l, char d, void *data);

static FLT_OR_DBL
expAptamerContrib(int i, int j, int k, int l, char d, void *data);

static FLT_OR_DBL
expAptamerContribHairpin(int i, int j, int k, int l, char d, void *data);

static vrna_basepair_t *
backtrack_int_motif(int i, int j, int k, int l, char d, void *data);

static vrna_basepair_t *
backtrack_hp_motif(int i, int j, int k, int l, char d, void *data);

static quadruple_position *
scanForMotif( const char *seq,
              const char *motif1,
              const char *motif2);

static vrna_basepair_t *
scanForPairs( const char  *motif5,
              const char  *motif3,
              int         *pair_count);

static int
vrna_sc_detect_hi_motif(vrna_fold_compound_t *vc,
                        const char *structure,
                        int *i,
                        int *j,
                        int *k,
                        int *l);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC vrna_sc_motif_t *
vrna_sc_ligand_detect_motifs( vrna_fold_compound_t *vc,
                              const char *structure){

  int             a, b, c, d, motif_count, motif_list_size;
  vrna_sc_motif_t *motif_list;
  
  motif_list = NULL;
  
  if(vc && structure){
    a               = 1;
    motif_count     = 0;
    motif_list_size = 10;

    motif_list = (vrna_sc_motif_t *)vrna_alloc(sizeof(vrna_sc_motif_t) * motif_list_size);

    while(vrna_sc_detect_hi_motif(vc, structure, &a, &b, &c, &d)){
      if(motif_count == motif_list_size){
        motif_list_size *= 1.2;
        motif_list = (vrna_sc_motif_t *)vrna_realloc(motif_list, sizeof(vrna_sc_motif_t) * motif_list_size);
      }

      motif_list[motif_count].number  = 0;

      if(c != 0){ /* interior loop motif */
        motif_list[motif_count].i = a;
        motif_list[motif_count].j = b;
        motif_list[motif_count].k = c;
        motif_list[motif_count].l = d;
        a = c;
      } else { /* hairpin loop motif */
        motif_list[motif_count].i = a;
        motif_list[motif_count].j = b;
        motif_list[motif_count].k = a;
        motif_list[motif_count].l = b;
        a = b;
      }
      motif_count++;
    }

    motif_list = (vrna_sc_motif_t *)vrna_realloc(motif_list, sizeof(vrna_sc_motif_t) * (motif_count + 1));
    motif_list[motif_count].i       = 0;
    motif_list[motif_count].j       = 0;
    motif_list[motif_count].k       = 0;
    motif_list[motif_count].l       = 0;
  }

  return motif_list;
}

PRIVATE int
vrna_sc_detect_hi_motif(vrna_fold_compound_t *vc,
                        const char *structure,
                        int *i,
                        int *j,
                        int *k,
                        int *l){

  int p, q, n;
  quadruple_position  *pos;
  ligand_data         *ldata;

  if(vc && vc->sc && vc->sc->data){
    n = vc->length;
    ldata = (ligand_data *)vc->sc->data;

    for(p = *i; p < n; p++){
      for(pos = ldata->positions; pos->i; pos++){
        if(pos->i == p){
          /* check whether we find the motif in the provided structure */
          int i_m, j_m, k_m, l_m;
          i_m = pos->i;
          j_m = pos->j;
          k_m = pos->k;
          l_m = pos->l;
          for(q = 0; q < strlen(ldata->struct_motif_5); q++){
            if(ldata->struct_motif_5[q] != structure[i_m+q-1])
              break;
          }
          if(q == strlen(ldata->struct_motif_5)){ /* 5' motif detected */
            if(k_m > 0){
              for(q = 0; q < strlen(ldata->struct_motif_3); q++){
                if(ldata->struct_motif_3[q] != structure[l_m+q-1])
                  break;
              }
              if(q == strlen(ldata->struct_motif_3)){ /* 3' motif detected */
                *i = i_m;
                *j = j_m;
                *k = k_m;
                *l = l_m;
                return 1;
              }
            } else {
                *i = i_m;
                *j = j_m;
                *k = k_m;
                *l = l_m;
                return 1;
            }
          }
        }
      }
    }
      
  }
  return 0;
}

PUBLIC int
vrna_sc_get_hi_motif( vrna_fold_compound_t *vc,
                      int *i,
                      int *j,
                      int *k,
                      int *l){

  int p, n;
  quadruple_position  *pos;
  ligand_data         *ldata;

  if(vc && vc->sc && vc->sc->data){
    n = vc->length;
    ldata = (ligand_data *)vc->sc->data;

    for(p = *i; p < n; p++){
      for(pos = ldata->positions; pos->i; pos++){
        if(pos->i == p){
          *i = pos->i;
          *j = pos->j;
          *k = pos->k;
          *l = pos->l;
          return 1;
        }
      }
    }
  }
  return 0;
}

PUBLIC int
vrna_sc_add_hi_motif( vrna_fold_compound_t *vc,
                      const char *seq,
                      const char *structure,
                      FLT_OR_DBL energy,
                      unsigned int options){

    int                   i, cp, cp2;
    char                  *sequence, *motif, *motif_alt;
    vrna_md_t             *md_p;
    ligand_data           *ldata;

    sequence              = NULL;
    motif                 = NULL;
    motif_alt             = NULL;
    ldata                 = NULL;
    md_p                  = NULL;

    sequence  = vrna_cut_point_remove(seq, &cp);                /* ligand sequence motif  */
    motif     = vrna_cut_point_remove(structure, &cp2);         /* ligand structure motif */

    /* check for obvious inconsistencies in input sequence/structure motif */
    if(cp != cp2){
      vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: Cutpoint in sequence and structure motif differ!");
      goto hi_motif_error;
    } else if(strlen(seq) != strlen(structure)){
      vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: length of sequence and structure motif differ!");
      goto hi_motif_error;
    }

    /* create auxiliary soft constraints data structure */
    ldata                 = vrna_alloc(sizeof(ligand_data));
    ldata->seq_motif_5    = NULL;
    ldata->seq_motif_3    = NULL;
    ldata->struct_motif_5 = NULL;
    ldata->struct_motif_3 = NULL;
    ldata->positions      = NULL;
    ldata->energy         = (int)(energy * 100.);

    split_sequence(sequence, &(ldata->seq_motif_5), &(ldata->seq_motif_3), cp);
    split_sequence(motif, &(ldata->struct_motif_5), &(ldata->struct_motif_3), cp);

    motif_alt = vrna_alloc(sizeof(char) * (strlen(motif) + 1)); /* alternative structure motif */
    memset(motif_alt, '.', strlen(motif) - 1);

    if(cp > 0){
      if((motif[0] != '(') || (motif[strlen(motif) - 1] != ')') || (motif[cp-2] != '(') || (motif[cp-1] != ')')){
        vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: No closing and/or enclosed pair in interior loop motif!");
        goto hi_motif_error;
      }
      /* construct corresponding alternative interior loop motif (....(&)...) */
      motif_alt[0] = '(';
      motif_alt[cp-2] = '(';
      motif_alt[cp-1] = ')';
      motif_alt[strlen(motif) - 1] = ')';
      motif_alt[strlen(motif)] = '\0';

      vrna_sc_add_bt(vc, &backtrack_int_motif);
      vrna_sc_add_f(vc, &AptamerContrib);
      vrna_sc_add_exp_f(vc, &expAptamerContrib);

    } else {
      if((motif[0] != '(') || (motif[strlen(motif) - 1] != ')')){
        vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: No closing pair in hairpin motif!");
        goto hi_motif_error;
      }

      /* construct corresponding alternative hairpin motif (....) */
      motif_alt[0] = '(';
      motif_alt[strlen(motif) - 1] = ')';
      motif_alt[strlen(motif)] = '\0';

      vrna_sc_add_bt(vc, &backtrack_hp_motif);
      vrna_sc_add_f(vc, &AptamerContribHairpin);
      vrna_sc_add_exp_f(vc, &expAptamerContribHairpin);
    }

    /* correct motif contributions */
    if(vc->params)
      md_p = &(vc->params->model_details);
    else
      md_p = &(vc->exp_params->model_details);

    correctMotifContribution(seq, motif, motif_alt, &(ldata->energy), &(ldata->energy_alt), md_p);

    /* scan for sequence motif positions */
    ldata->positions = scanForMotif(vc->sequence, ldata->seq_motif_5, ldata->seq_motif_3);

    /* scan for additional base pairs in the structure motif */
    int pair_count = 0;
    vrna_basepair_t *pairs = scanForPairs(ldata->struct_motif_5, ldata->struct_motif_3, &pair_count);
    if((pair_count > 0) && (pairs == NULL)){ /* error while parsing structure motif */
      vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: Error while parsing additional pairs in structure motif");
      goto hi_motif_error;
    }

    ldata->pairs      = pairs;
    ldata->pair_count = pair_count;

    /* add generalized soft-constraint data structure and corresponding 'delete' function */
    vrna_sc_add_data(vc, (void *)ldata, &delete_ligand_data);

    free(sequence);
    free(motif);
    free(motif_alt);

    return 1; /* success */

/* exit with error */
hi_motif_error:

    free(sequence);
    free(motif);
    free(motif_alt);
    delete_ligand_data(ldata);

    return 0;
}

static void
split_sequence( const char *string,
                      char **seq1,
                      char **seq2,
                      int cp){

  int l = (int)strlen(string);
  *seq1 = NULL;
  *seq2 = NULL;

  if(cp > 0){
    if(cp < l){
      *seq1 = vrna_alloc(sizeof(char) * cp);
      strncpy(*seq1, string, cp - 1);
      (*seq1)[cp - 1] = '\0';
      *seq2 = vrna_alloc(sizeof(char) * (l - cp + 2));
      strncpy(*seq2, string + cp - 1, (l - cp + 1));
      (*seq2)[l - cp + 1] = '\0';
    }
  } else {
    *seq1 = vrna_alloc(sizeof(char) * (l+1));
    strncpy(*seq1, string, l);
    (*seq1)[l] = '\0';
  }
}

static void
correctMotifContribution( const char *seq,
                          const char *struct_motif,
                          const char *struct_motif_alt,
                          int *contribution,
                          int *contribution_alt,
                          vrna_md_t *md){

  float                 alt, corr, energy;
  vrna_fold_compound_t  *tmp_vc;

  tmp_vc  = vrna_fold_compound(seq, md, VRNA_OPTION_EVAL_ONLY);
  alt     = vrna_eval_structure(tmp_vc, struct_motif_alt);
  corr    = vrna_eval_structure(tmp_vc, struct_motif);
  energy  = corr - alt;

  *contribution     += (int)(energy * 100.);
  *contribution_alt  = (int)(alt    * 100.);

  vrna_fold_compound_free(tmp_vc);
}

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

static int
AptamerContrib(int i, int j, int k, int l, char d, void *data){

  quadruple_position  *pos;
  ligand_data         *ldata;

  if(d == VRNA_DECOMP_PAIR_IL){
    ldata = (ligand_data *)data;
    for(pos = ((ligand_data *)data)->positions; pos->i; pos++){
      if((pos->i == i) && (pos->j == j) && (pos->k == k) && (pos->l == l)){
        return ldata->energy;
      }
    }
  }

  return 0;
}

static int
AptamerContribHairpin(int i, int j, int k, int l, char d, void *data){

  quadruple_position  *pos;
  ligand_data         *ldata;

  if(d == VRNA_DECOMP_PAIR_HP){
    ldata = (ligand_data *)data;
    for(pos = ((ligand_data *)data)->positions; pos->i; pos++){
      if((pos->i == i) && (pos->j == j)){
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

  exp_e = 1.;

  if(d == VRNA_DECOMP_PAIR_IL){
    ldata = (ligand_data *)data;
    kT    = (37. + K0) * GASCONST;

    for(pos = ldata->positions; pos->i; pos++){
      if((pos->i == i) && (pos->j == j) && (pos->k == k) && (pos->l == l)){
        exp_e =   (FLT_OR_DBL)exp((double) (-ldata->energy) * 10./kT);
        exp_e +=  (FLT_OR_DBL)exp((double) (-ldata->energy_alt) * 10./kT); /* add alternative, i.e. unbound ligand */
        break;
      }
    }
  }

  return exp_e;
}

static FLT_OR_DBL
expAptamerContribHairpin(int i, int j, int k, int l, char d, void *data){

  quadruple_position  *pos;
  ligand_data         *ldata;
  FLT_OR_DBL          exp_e;
  double              kT;

  exp_e = 1.;

  if(d == VRNA_DECOMP_PAIR_HP){
    ldata = (ligand_data *)data;
    kT    = (37. + K0) * GASCONST;

    for(pos = ldata->positions; pos->i; pos++){
      if((pos->i == i) && (pos->j == j)){
        exp_e =   (FLT_OR_DBL)exp((double) (-ldata->energy) * 10./kT);
        exp_e +=  (FLT_OR_DBL)exp((double) (-ldata->energy_alt) * 10./kT); /* add alternative, i.e. unbound ligand */
        break;
      }
    }
  }

  return exp_e;
}

static vrna_basepair_t *
backtrack_int_motif(int i, int j, int k, int l, char d, void *data){

  int                 bp_size = 15;
  vrna_basepair_t     *pairs = NULL;
  quadruple_position  *pos;
  ligand_data         *ldata;

  if(d == VRNA_DECOMP_PAIR_IL){
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
        for(count = 0,pptr = ldata->pairs; pptr && (pptr->i != 0); pptr++, count++){
          pairs[count].i = (pptr->i < 0) ? j + pptr->i : i + pptr->i - 1;
          pairs[count].j = (pptr->j < 0) ? j + pptr->j : i + pptr->j - 1;
        }
        pairs[count].i = pairs[count].j = 0;
#endif

        return pairs;
      }
    }
  }

  return pairs;
}

static vrna_basepair_t *
backtrack_hp_motif(int i, int j, int k, int l, char d, void *data){

  int                 count;
  vrna_basepair_t     *pairs = NULL;
  quadruple_position  *pos;
  ligand_data         *ldata;
  vrna_basepair_t     *pptr;

  if(d == VRNA_DECOMP_PAIR_HP){
    ldata = (ligand_data *)data;
    for(pos = ldata->positions; pos->i; pos++){
      if((pos->i == i) && (pos->j == j)){
        /* found motif in our list, lets create pairs */
        pairs = vrna_alloc(sizeof(vrna_basepair_t) * (ldata->pair_count + 1));
        for(count = 0,pptr = ldata->pairs; pptr && (pptr->i != 0); pptr++, count++){
          pairs[count].i = i + pptr->i - 1;
          pairs[count].j = i + pptr->j - 1;
        }
        pairs[count].i = pairs[count].j = 0;
        return pairs;
      }
    }
  }

  return pairs;
}

static quadruple_position *
scanForMotif( const char *seq,
              const char *motif1,
              const char *motif2){

  int   i, j, k, l, l1, l2, n, cnt, cnt2;
  char  *ptr;
  quadruple_position *pos;
  
  n     = (int) strlen(seq);
  l1    = (int) strlen(motif1);
  l2    = (motif2) ? (int) strlen(motif2) : 0;
  cnt   = 0;
  cnt2  = 5; /* initial guess how many matching motifs we might encounter */

  pos = (quadruple_position *)vrna_alloc(sizeof(quadruple_position) * cnt2);

  for(i = 0; i <= n - l1 - l2; i++){
    if(seq[i] == motif1[0]){
      for(j = i+1; j < i + l1; j++){
        if(seq[j] == motif1[j-i]){
          continue;
        }
        else goto next_i;
      }
      /* found 5' motif */
      if(motif2){
        for(k = j + 1; k <= n - l2; k++){
          if(seq[k] == motif2[0]){
            for(l = k + 1; l < k + l2; l++){
              if(seq[l] == motif2[l-k]){
                continue;
              }
              else goto next_k;
            }
            /* we found a quadruple, so store it */
            pos[cnt].i   = i + 1;
            pos[cnt].j   = l;
            pos[cnt].k   = j;
            pos[cnt++].l = k + 1;

            /* allocate more memory if necessary */
            if(cnt == cnt2){
              cnt2 *= 2;
              pos = (quadruple_position *)vrna_realloc(pos, sizeof(quadruple_position) * cnt2);
            }
          }
/* early exit from l loop */
next_k: continue;
        }
      } else { /* hairpin loop motif */
        /* store it */
        pos[cnt].i   = i + 1;
        pos[cnt].j   = j;
        pos[cnt].k   = 0;
        pos[cnt++].l = 0;

        /* allocate more memory if necessary */
        if(cnt == cnt2){
          cnt2 *= 2;
          pos = (quadruple_position *)vrna_realloc(pos, sizeof(quadruple_position) * cnt2);
        }
      }
    }
/* early exit from j loop */
next_i: continue;
  }

  /* reallocate to actual size */
  pos = (quadruple_position *)vrna_realloc(pos, sizeof(quadruple_position) * (cnt + 1));

  /* set end marker */
  pos[cnt].i = pos[cnt].j = pos[cnt].k = pos[cnt].l = 0;

  return pos;
}

static vrna_basepair_t *
scanForPairs( const char  *motif5,
              const char  *motif3,
              int         *pair_count){

  int             i, l5, l3, stack_size, stack_count, *stack;
  vrna_basepair_t *pairs;

  l5          = (motif5) ? strlen(motif5) : 0;
  l3          = (motif3) ? strlen(motif3) : 0;
  stack_count = 0;
  stack_size  = l5 + l3 + 1;
  *pair_count = 0;
  stack       = vrna_alloc(sizeof(int)              * stack_size);
  pairs       = vrna_alloc(sizeof(vrna_basepair_t)  * stack_size);

  /* go through 5' side of structure motif */
  for(i = 2; i < l5; i++){
    if(motif5[i - 1] == '('){
      stack[stack_count++] = i;
    } else if(motif5[i - 1] == ')'){
      pairs[*pair_count].i = stack[--stack_count];
      pairs[*pair_count].j = i;
      /* printf("5' p[%d, %d]\n", pairs[*pair_count].i, pairs[*pair_count].j); */
      (*pair_count)++;
      if(stack_count < 0){
        vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: 5' structure motif contains unbalanced brackets");
        free(stack);
        free(pairs);
        return NULL;
      }
    }
  }

  if(motif3){
    for(i = 2; i < l3; i++){ /* go through 3' side of motif */
      if(motif3[i-1] == '('){
        stack[stack_count++] = -(l3 - i);
      } else if(motif3[i-1] == ')'){
        pairs[*pair_count].i = stack[--stack_count];
        pairs[*pair_count].j = -(l3 - i);
        /* printf("3' p[%d, %d]\n", pairs[*pair_count].i, pairs[*pair_count].j); */
        (*pair_count)++;
        if(stack_count < 0){
          vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: 3' structure motif contains unbalanced brackets");
          free(stack);
          free(pairs);
          return NULL;
        }
      }
    }
  }

  if(stack_count != 0){
    vrna_message_warning("vrna_sc_add_ligand_binding@ligand.c: structure motif contains unbalanced brackets");
    (*pair_count)++;
    free(stack);
    free(pairs);
    return NULL;
  }

  if(*pair_count > 0){
    pairs = vrna_realloc(pairs, sizeof(vrna_basepair_t) * (*pair_count + 1));
    pairs[*pair_count].i = pairs[*pair_count].j = 0;
  } else {
    free(pairs);
    pairs = NULL;
  }

  free(stack);

  return pairs;
}

