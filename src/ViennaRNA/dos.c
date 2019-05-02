#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/mfe.h"
#include "ViennaRNA/file_utils.h"
#include "ViennaRNA/datastructures/hash_tables.h"

struct aux_arrays {
  int *cc; /* auxilary arrays for canonical structures     */
  int *cc1; /* auxilary arrays for canonical structures     */
  int *Fmi; /* holds row i of fML (avoids jumps in memory)  */
  int *DMLi; /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  int *DMLi1; /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  int *DMLi2; /*                MIN(fML[i+2,k]+fML[k+1,j])    */
};

PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int length)
{
  unsigned int j;
  struct aux_arrays *aux = (struct aux_arrays *) vrna_alloc (sizeof(struct aux_arrays));

  aux->cc = (int *) vrna_alloc (sizeof(int) * (length + 2)); /* auxilary arrays for canonical structures     */
  aux->cc1 = (int *) vrna_alloc (sizeof(int) * (length + 2)); /* auxilary arrays for canonical structures     */
  aux->Fmi = (int *) vrna_alloc (sizeof(int) * (length + 1)); /* holds row i of fML (avoids jumps in memory)  */
  aux->DMLi = (int *) vrna_alloc (sizeof(int) * (length + 1)); /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  aux->DMLi1 = (int *) vrna_alloc (sizeof(int) * (length + 1)); /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  aux->DMLi2 = (int *) vrna_alloc (sizeof(int) * (length + 1)); /*                MIN(fML[i+2,k]+fML[k+1,j])    */

  /* prefill helper arrays */
  for (j = 0; j <= length; j++)
    aux->Fmi[j] = aux->DMLi[j] = aux->DMLi1[j] = aux->DMLi2[j] = INF;

  return aux;
}

PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux, unsigned int length)
{
  unsigned int j;
  int *FF;

  FF = aux->DMLi2;
  aux->DMLi2 = aux->DMLi1;
  aux->DMLi1 = aux->DMLi;
  aux->DMLi = FF;
  FF = aux->cc1;
  aux->cc1 = aux->cc;
  aux->cc = FF;
  for (j = 1; j <= length; j++)
    aux->cc[j] = aux->Fmi[j] = aux->DMLi[j] = INF;
}

PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux)
{
  free (aux->cc);
  free (aux->cc1);
  free (aux->Fmi);
  free (aux->DMLi);
  free (aux->DMLi1);
  free (aux->DMLi2);
  free (aux);
}

typedef struct key_value_ {
  unsigned long key;
  unsigned long value;
} key_value;

static unsigned
hash_function_dos(void *hash_entry, unsigned long hashtable_size)
{
  key_value *hem = ((key_value *) hash_entry);
  unsigned long c = (unsigned long)hem->key;
  return c % hashtable_size;
}

/***
 * 0 is equal, 1 is different!
 */
static int
hash_comparison_dos(void *x, void *y)
{
  key_value *hem_x = ((key_value *) x);
  key_value *hem_y = ((key_value *) y);

  if ((x == NULL) ^ (y == NULL))
    return 1;

  return !(hem_x->key == hem_y->key);
}

static int
free_hash_entry_dos(void *hash_entry)
{
  /*
  if (hash_entry != NULL) {
    key_value *hem_x = ((key_value *) hash_entry);
   free (hem_x);
  }
  */
  return 0;
}

static vrna_hash_table_t
create_hashtable()
{
  vrna_callback_ht_free_entry *my_free = free_hash_entry_dos;
  vrna_callback_ht_compare_entries *my_comparison = hash_comparison_dos;
  vrna_callback_ht_hash_function *my_hash_function = hash_function_dos;
  vrna_hash_table_t ht = vrna_ht_init (20, my_comparison, my_hash_function, my_free);
  return ht;
}

typedef struct energy_count_ {
  int energy;
  unsigned long count;
} energy_count;

typedef struct hashtable_list_ {
  unsigned long length;
  unsigned long allocated_size;
  energy_count * list_energy_count_pairs;
  key_value ** list_key_value_pairs;
  vrna_hash_table_t ht_energy_index; // lookup table;
} hashtable_list;

struct dp_counts_per_energy {
  hashtable_list *n_ij_e;
  hashtable_list *n_ij_A_e;
  hashtable_list *n_ij_M_e;
  hashtable_list *n_ij_M1_e;
};

static hashtable_list
create_hashtable_list()
{
  hashtable_list ht_list;
  ht_list.allocated_size = 10;
  ht_list.length = 0;
  ht_list.list_energy_count_pairs = vrna_alloc (sizeof(energy_count) * ht_list.allocated_size);
  ht_list.list_key_value_pairs = vrna_alloc (sizeof(key_value*) * ht_list.allocated_size);
  ht_list.ht_energy_index = create_hashtable ();
  return ht_list;
}

static
void free_hashtable_list(hashtable_list *ht_list)
{
  vrna_ht_free(ht_list->ht_energy_index);
  free(ht_list->list_energy_count_pairs);

  int i=0;
  for(; i < ht_list->length; i++)
    free(ht_list->list_key_value_pairs[i]);

  free(ht_list->list_key_value_pairs);
}

static
void
hashtable_list_add_count(hashtable_list *htl, int energy, int count)
{
  if (htl->ht_energy_index != NULL) {
    key_value to_check;
    to_check.key = energy;
    //to_check.value = count;
    key_value* lookup_result = NULL;
    lookup_result = vrna_ht_get (htl->ht_energy_index, (void *) &to_check);
    if (lookup_result == NULL) {
      //value is not in list.
      if (htl->length >= (htl->allocated_size-1)) {
        htl->allocated_size += 10;
        htl->list_energy_count_pairs = vrna_realloc (htl->list_energy_count_pairs, sizeof(energy_count) * htl->allocated_size);
        htl->list_key_value_pairs = vrna_realloc (htl->list_key_value_pairs, sizeof(key_value*) * htl->allocated_size);
      }
      energy_count ec;
      ec.count = count;
      ec.energy = energy;
      unsigned long list_index = htl->length;
      htl->list_energy_count_pairs[list_index] = ec;
      to_check.value = list_index;
      key_value *to_store = vrna_alloc(sizeof(key_value));
      *to_store = to_check;
      htl->list_key_value_pairs[list_index] = to_store;
      htl->length++;

      key_value* to_insert = htl->list_key_value_pairs[list_index];
      int res = vrna_ht_insert (htl->ht_energy_index, (void *) to_insert);
      if (res != 0) {
        fprintf (stderr, "dos.c: hash table insert failed!");
      }
    }
    else {
      // the energy-index pair is already in the list.
      unsigned long list_index = lookup_result->value;
      htl->list_energy_count_pairs[list_index].count += count;
    }
  }
}

static
unsigned long
lookup_count_from_energy(hashtable_list *htl, int energy)
{
  key_value to_check;
  to_check.key = energy;
  key_value* lookup_result;
  lookup_result = vrna_ht_get (htl->ht_energy_index, (void *) &to_check);
  if (lookup_result == NULL) {
    //value is not in list.
    return 0;
  }
  else {
    unsigned long index = lookup_result->value;
    return htl->list_energy_count_pairs[index].count;
  }
}

static
unsigned long
lookup_count_from_index(hashtable_list *htl, unsigned long index)
{
  if (index < 0)
    return 0;
  return htl->list_energy_count_pairs[index].count;
}

int
lookup_energy_from_index(hashtable_list *htl, unsigned long index)
{
  if (index < 0)
    return 1000000;
  return htl->list_energy_count_pairs[index].energy;
}

PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc, int i, int j, struct aux_arrays *aux, int min_energy, int max_e, int step_energy, int energy_length,
               struct dp_counts_per_energy *dp_count_matrix_pt)
{
  unsigned char hc_decompose;
  unsigned int n;
  int e, new_c, energy, stackEnergy, ij, dangle_model, noLP, *DMLi1, *DMLi2, *cc, *cc1;
  hashtable_list *source_table;
  hashtable_list *source_table_2;
  hashtable_list *result_table;

  n = fc->length;
  ij = fc->jindx[j] + i;
  dangle_model = fc->params->model_details.dangles;
  noLP = fc->params->model_details.noLP;
  hc_decompose = fc->hc->mx[n * i + j];
  DMLi1 = aux->DMLi1;
  DMLi2 = aux->DMLi2;
  cc = aux->cc;
  cc1 = aux->cc1;
  e = INF;
  int *rtype = &(fc->params->model_details.rtype[0]);
  short *S1 = fc->sequence_encoding;

  int type = fc->ptype[fc->jindx[j] + i];

  int no_close = (((type == 3) || (type == 4)) && fc->params->model_details.noGUclosure);

  /* do we evaluate this pair? */
  if (type) { //(hc_decompose) {
    new_c = INF;

    /* check for hairpin loop */
    int energy_hp = E_Hairpin (j - i - 1, type, S1[i + 1], S1[j - 1], fc->sequence + i - 1, fc->params); // vrna_E_hp_loop(fc, i, j);
    fc->matrices->c[ij] = MIN2(fc->matrices->c[ij], energy_hp);
    new_c = MIN2(new_c, energy_hp);

    int max_energy = min_energy + (energy_length - 1) * step_energy;

    if (energy_hp <= max_e) {
      //dp_count_matrix_pt->n_ij_e[ij][energy_hp]++;
      result_table = &dp_count_matrix_pt->n_ij_e[ij];
      hashtable_list_add_count (result_table, energy_hp, 1);
    }

    /* check for interior loops */
    int min_int = new_c;
    int maxp = MIN2(j - 2 - TURN, i + MAXLOOP + 1);
    int p, q, type_2, pq;
    for (p = i + 1; p <= maxp; p++) {
      unsigned int minq = p + TURN + 1;

      for (q = minq; q < j; q++) {
        pq = fc->jindx[q] + p;
        /* set distance to reference structure... */
        type_2 = fc->ptype[fc->jindx[q] + p];

        if (type_2 == 0)
          continue;

        type_2 = rtype[type_2];

        if (no_closingGU)
          if (no_close || (type_2 == 3) || (type_2 == 4))
            if ((p > i + 1) || (q < j - 1))
              continue;

        /* continue unless stack */
        energy = E_IntLoop (p - i - 1, j - q - 1, type, type_2, S1[i + 1], S1[j - 1], S1[p - 1], S1[q + 1], fc->params);
        min_int = MIN2(min_int, fc->matrices->c[pq] + energy);
        if (fc->matrices->c[pq] != INF) {
          fc->matrices->c[ij] = MIN2(fc->matrices->c[ij], fc->matrices->c[pq] + energy);

          source_table = &dp_count_matrix_pt->n_ij_e[pq];
          result_table = &dp_count_matrix_pt->n_ij_e[ij];
          for (int ei = 0; ei < source_table->length; ei++) {
            int pq_energy = lookup_energy_from_index (source_table, ei);
            unsigned long pq_count = lookup_count_from_index (source_table, ei);
            int sum_energy = pq_energy + energy;
            if ((sum_energy >= min_energy) && (sum_energy <= max_e))
              hashtable_list_add_count (result_table, sum_energy, pq_count);
          }
        }
      } /* end q-loop */
    } /* end p-loop */
    new_c = MIN2(new_c, min_int);

    /* check for multibranch loops */
    int min_ml = INF;
    if (!no_close) {
      /* dangle energies for multiloop closing stem */
      int tt = rtype[type];
      int temp2 = fc->params->MLclosing;
      if (dangles == 2)
        temp2 += E_MLstem (tt, S1[j - 1], S1[i + 1], fc->params);
      else
        temp2 += E_MLstem (tt, -1, -1, fc->params);

      int u;
      for (u = i + TURN + 2; u < j - TURN - 2; u++) {
        int i1u = fc->jindx[u] + (i + 1);
        int u1j1 = fc->jindx[j - 1] + (u + 1);

        if ((fc->matrices->fML[i1u] != INF) && (fc->matrices->fM1[u1j1] != INF)) {
          int tmp_min_ml = MIN2(fc->matrices->c[ij], fc->matrices->fML[i1u] + fc->matrices->fM1[u1j1] + temp2);
          min_ml = MIN2(min_ml, tmp_min_ml);

          fc->matrices->c[ij] = MIN2(fc->matrices->c[ij], fc->matrices->fML[i1u] + fc->matrices->fM1[u1j1] + temp2);

          source_table = &dp_count_matrix_pt->n_ij_M_e[i1u];
          source_table_2 = &dp_count_matrix_pt->n_ij_M1_e[u1j1];
          result_table = &dp_count_matrix_pt->n_ij_e[ij];
          for (int ei_1 = 0; ei_1 < source_table->length; ei_1++) {
            int m_energy = lookup_energy_from_index (source_table, ei_1);
            unsigned long m_count = lookup_count_from_index (source_table, ei_1);
            for (int ei_2 = 0; ei_2 < source_table_2->length; ei_2++) {
              int m1_energy = lookup_energy_from_index (source_table_2, ei_2);
              unsigned long m1_count = lookup_count_from_index (source_table_2, ei_2);
              int sum_energy = m_energy + m1_energy + temp2;
              if ((sum_energy >= min_energy) && (sum_energy <= max_e)) {
                unsigned long c_count = m_count * m1_count;
                hashtable_list_add_count(result_table, sum_energy, c_count);
              }
            }
          }
        }
      }
    }

    new_c = MIN2(new_c, min_ml);

    if (dangle_model == 3) {
      /* coaxial stacking */
      energy = vrna_E_mb_loop_stack (fc, i, j);
      new_c = MIN2(new_c, energy);
    }

    /* remember stack energy for --noLP option */
    if (noLP) {
      stackEnergy = vrna_E_stack (fc, i, j);
      new_c = MIN2(new_c, cc1[j - 1] + stackEnergy);
      cc[j] = new_c;
      if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
        cc[j] -= fc->pscore[ij];

      e = cc1[j - 1] + stackEnergy;
    }
    else {
      e = new_c;
    }

    /* finally, check for auxiliary grammar rule(s) */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_c)) {
      energy = fc->aux_grammar->cb_aux_c (fc, i, j, fc->aux_grammar->data);
      new_c = MIN2(new_c, energy);
    }

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      e -= fc->pscore[ij];
  } /* end >> if (pair) << */

  return e;
}

int
print_matrix(int *matrix, int *jiindx, int length_row, int length_col)
{
  int i, j;
  printf ("\n");
  for (i = 1; i <= length_row; i++) {
    for (j = 1; j <= length_col; j++) {
      int mj = fmax (i, j);
      int mi = fmin (i, j);
      int ij = jiindx[mj] + mi;
      printf (" %7.4g", (float) matrix[ij]);
    }
    printf ("\n");
  }

  return 0;
}

int
print_matrix_energy_counts(int **matrix, int *jiindx, int length_row, int length_col, int min_energy, int step_energy, int energy_length)
{
  int i, j;
  printf ("\n");
  for (i = 1; i <= length_row; i++) {
    for (j = i; j <= length_col; j++) {
      int mj = fmax (i, j);
      int mi = fmin (i, j);
      int ij = jiindx[mj] + mi;
      if (matrix[ij] != NULL) {
        int e_idx = 0;
        for (; e_idx < energy_length; e_idx++) {
          unsigned long e_count = matrix[ij][e_idx];
          int e = min_energy + e_idx * step_energy;
          if (e_count > 0)
            printf (" %d %d %.2f %d\n", i, j, e / 100., e_count);
        }
      }
    }
  }
  printf ("\n");
  return 0;
}

int
print_array_energy_counts(hashtable_list *matrix, int length_col, int min_energy, int step_energy, int energy_length)
{
  int i, j;
  printf ("\n");
  int max_energy = min_energy + energy_length - 1;
  printf ("min_e: %d, max_e: %d, length: %d, energy_length: %d\n", min_energy, max_energy, length_col, energy_length);
  j = length_col;
  if (matrix[j].length > 0){
    for (int e = min_energy; e <= max_energy; e++){
      unsigned long count = lookup_count_from_energy(&matrix[j], e);
      if (count > 0){
        printf ("%6.2f\t%d\n", e / 100., count);
      }
    }
  }
  printf ("\n");
  return 0;
}

/* fill DP matrices */
PRIVATE int
fill_arrays(vrna_fold_compound_t *fc)
{
  int i, j, ij, length, turn, uniq_ML, *indx, *f5, *c, *fML, *fM1;
  vrna_param_t *P;
  vrna_mx_mfe_t *matrices;
  vrna_ud_t *domains_up;
  struct aux_arrays *helper_arrays;

  length = (int) fc->length;
  indx = fc->jindx;
  P = fc->params;
  uniq_ML = P->model_details.uniq_ML;
  turn = P->model_details.min_loop_size;
  matrices = fc->matrices;
  f5 = matrices->f5;
  c = matrices->c;
  fML = matrices->fML;
  fM1 = matrices->fM1;
  domains_up = fc->domains_up;
  int dangles = P->model_details.dangles;
  short *S1 = fc->sequence_encoding;
  int circ = P->model_details.circ;

  hashtable_list *source_table;
  hashtable_list *source_table_2;
  hashtable_list *result_table;

  hashtable_list *n_ij_e = (hashtable_list *) vrna_alloc (sizeof(hashtable_list) * ((length + 1) * (length + 1)));
  hashtable_list *n_ij_A_e = (hashtable_list *) vrna_alloc (sizeof(hashtable_list) * (length + 1));
  hashtable_list *n_ij_M_e = (hashtable_list *) vrna_alloc (sizeof(hashtable_list) * ((length + 1) * (length + 1)));
  hashtable_list *n_ij_M1_e = (hashtable_list *) vrna_alloc (sizeof(hashtable_list) * ((length + 1) * (length + 1)));

  struct dp_counts_per_energy count_matrix_pt;
  count_matrix_pt.n_ij_e = n_ij_e;
  count_matrix_pt.n_ij_A_e = n_ij_A_e;
  count_matrix_pt.n_ij_M_e = n_ij_M_e;
  count_matrix_pt.n_ij_M1_e = n_ij_M1_e;

  //printf("mfe %f\n", vrna_mfe(fc, NULL));

  int min_energy = (int) round (vrna_mfe (fc, NULL) * 100.0); // -500; // TODO compute mfe.
  printf ("min_energy (global): %d \n", min_energy);
  /* search through DP matrices for minimal entry */
  for (int i = 1; i < length; i++)
    for (int j = i + 1; j <= length; j++) {
      ij = indx[j] + i;
      int e = fc->matrices->c[ij];
      if (e < min_energy)
        min_energy = e;

      e = fc->matrices->fML[ij];
      if (e < min_energy)
        min_energy = e;
    }

  for (int i = 1; i <= length; i++) {
    int e = fc->matrices->f5[i];
    if (e < min_energy)
      min_energy = e;
  }

  //min_energy -= 1000;
  int max_energy = 6000; // TODO compute max_energy.
  int step_energy = 1; // 10 decakal. (smallest unit of energy computations)
  int range = max_energy - min_energy;
  int energy_length = ceil (range / (float) step_energy) + 1;
  printf ("min_energy: %d \n", min_energy);
  printf ("max_energy: %d %d\n", max_energy, min_energy + (step_energy * (energy_length - 1)));
  printf ("range: %d, energy_length: %d\n", range, energy_length);

  /* allocate memory for all helper arrays */
  helper_arrays = get_aux_arrays (length);

  /* pre-processing ligand binding production rule(s) */
  if (domains_up && domains_up->prod_cb)
    domains_up->prod_cb (fc, domains_up->data);

  /* prefill matrices with init contributions */
  for (j = 1; j <= length; j++)
    for (i = (j > turn ? (j - turn) : 1); i <= j; i++) {
      c[indx[j] + i] = fML[indx[j] + i] = INF;
      if (uniq_ML)
        fM1[indx[j] + i] = INF;
    }

  /* start recursion */
  if (length <= turn) {
    /* clean up memory */
    free_aux_arrays (helper_arrays);

    /* return free energy of unfolded chain */
    return 0;
  }

  for (i = length - turn - 1; i >= 1; i--) {
    for (j = i + turn + 1; j <= length; j++) {
      ij = indx[j] + i;

      //printf ("alloc ij: %d %d %d\n", i, j, energy_length);

      int hc_decompose = fc->hc->mx[length * i + j];
      int type = fc->ptype[fc->jindx[j] + i];

      //prepare matrices (add third dimension)
      count_matrix_pt.n_ij_e[ij] = create_hashtable_list ();
      count_matrix_pt.n_ij_M_e[ij] = create_hashtable_list ();
      count_matrix_pt.n_ij_M1_e[ij] = create_hashtable_list ();

      /* do some pointer magic to allow negative array indices if necessary */
      // count_matrix_pt.n_ij_e[ij] -= min_energy;
      // count_matrix_pt.n_ij_M_e[ij] -= min_energy;
      // count_matrix_pt.n_ij_M1_e[ij] -= min_energy;
      /* decompose subsegment [i, j] with pair (i, j) */
      c[ij] = decompose_pair (fc, i, j, helper_arrays, min_energy, max_energy, step_energy, energy_length, &count_matrix_pt);

      /* decompose subsegment [i, j] that is multibranch loop part with at least one branch */
      int temp2 = 0;
      if (dangles == 2)
        temp2 += E_MLstem (type, (i == 1) ? S1[length] : S1[i - 1], S1[j + 1], P);
      else
        temp2 += E_MLstem (type, -1, -1, P);

      /* now to the actual computations... */
      /* 1st E_M[ij] = E_M1[ij] = E_C[ij] + b */
      int energy_index;
      if (matrices->c[ij] != INF) {
        matrices->fML[ij] = matrices->fM1[ij] = matrices->c[ij] + temp2;
        source_table = &count_matrix_pt.n_ij_e[ij];
        result_table = &count_matrix_pt.n_ij_M1_e[ij];
        hashtable_list *result_table_2 = &count_matrix_pt.n_ij_M_e[ij];
        for (int ei = 0; ei < source_table->length; ei++) {
          int c_energy = lookup_energy_from_index (source_table, ei);
          unsigned long c_count = lookup_count_from_index (source_table, ei);
          int sum_energy = c_energy + temp2;
          if ((sum_energy >= min_energy) && (sum_energy <= max_energy)){
            hashtable_list_add_count(result_table, sum_energy, c_count);
            hashtable_list_add_count(result_table_2, sum_energy, c_count);
          }
        }
      }
      else {
        matrices->fML[ij] = matrices->fM1[ij] = INF;
      }

      /* 2rd E_M[ij] = MIN(E_M[ij], E_M[i,j-1] + c) */
      if (fc->matrices->fML[fc->jindx[j - 1] + i] != INF) {
        matrices->fML[ij] = MIN2(matrices->fML[ij], matrices->fML[fc->jindx[j - 1] + i] + P->MLbase);

        source_table = &count_matrix_pt.n_ij_M_e[fc->jindx[j - 1] + i];
        result_table = &count_matrix_pt.n_ij_M_e[ij];
        for (int ei = 0; ei < source_table->length; ei++) {
          int m_energy = lookup_energy_from_index (source_table, ei);
          unsigned long m_count = lookup_count_from_index (source_table, ei);
          int sum_energy = m_energy + P->MLbase;
          if ((sum_energy >= min_energy) && (sum_energy <= max_energy)){
            hashtable_list_add_count(result_table, sum_energy, m_count);
          }
        }
      }

      /* 3th E_M1[ij] = MIN(E_M1[ij], E_M1[i,j-1] + c) */
      if (fc->matrices->fM1[fc->jindx[j - 1] + i] != INF) {
        matrices->fM1[ij] = MIN2(matrices->fM1[ij], matrices->fM1[fc->jindx[j - 1] + i] + P->MLbase);

        source_table = &count_matrix_pt.n_ij_M1_e[fc->jindx[j - 1] + i];
        result_table = &count_matrix_pt.n_ij_M1_e[ij];
        for (int ei = 0; ei < source_table->length; ei++) {
          int m1_energy = lookup_energy_from_index (source_table, ei);
          unsigned long m1_count = lookup_count_from_index (source_table, ei);
          int sum_energy = m1_energy + P->MLbase;
          if ((sum_energy >= min_energy) && (sum_energy <= max_energy)){
            hashtable_list_add_count(result_table, sum_energy, m1_count);
          }
        }
      }

      if (j > TURN + 2) {
        int u;
        int temp3;

        for (u = i; u < j; u++) {
          int u1j = fc->jindx[j] + u + 1;
          int iu = fc->jindx[u] + i;

          if (fc->matrices->c[u1j] != INF) {
            type = fc->ptype[u1j];

            /* [i..u] is unpaired */
            if (dangles == 2)
              temp2 = E_MLstem (type, S1[u], S1[j + 1], P);
            else
              temp2 = E_MLstem (type, -1, -1, P);

            temp3 = temp2 + (u - i + 1) * P->MLbase;
            matrices->fML[ij] = MIN2(matrices->fML[ij], matrices->c[u1j] + temp3);

        source_table = &count_matrix_pt.n_ij_e[u1j];
        result_table = &count_matrix_pt.n_ij_M_e[ij];
        for (int ei = 0; ei < source_table->length; ei++) {
          int c_energy = lookup_energy_from_index (source_table, ei);
          unsigned long c_count = lookup_count_from_index (source_table, ei);
          int sum_energy = c_energy + temp3;
          if ((sum_energy >= min_energy) && (sum_energy <= max_energy)){
            hashtable_list_add_count(result_table, sum_energy, c_count);
          }
        }

            /* [i...u] has at least one stem */
            if (fc->matrices->fML[iu] != INF) {
              matrices->fML[ij] = MIN2(matrices->fML[ij], matrices->fML[iu] + matrices->c[u1j] + temp2);

        source_table = &count_matrix_pt.n_ij_e[u1j];
        source_table_2 = &count_matrix_pt.n_ij_M_e[iu];
        result_table = &count_matrix_pt.n_ij_M_e[ij];
        for (int ei_1 = 0; ei_1 < source_table->length; ei_1++) {
          int c_energy = lookup_energy_from_index (source_table, ei_1);
          unsigned long c_count = lookup_count_from_index (source_table, ei_1);
          for (int ei_2 = 0; ei_2 < source_table_2->length; ei_2++) {
              int m_energy = lookup_energy_from_index (source_table_2, ei_2);
              unsigned long m_count = lookup_count_from_index (source_table_2, ei_2);
              int sum_energy = c_energy + m_energy + temp2;
              if ((sum_energy >= min_energy) && (sum_energy <= max_energy)){
                unsigned long product_count = c_count * m_count;
                hashtable_list_add_count(result_table, sum_energy, product_count);
              }
          }
        }
            }
          }
        }
      }

      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux))
        fc->aux_grammar->cb_aux (fc, i, j, fc->aux_grammar->data);
    } /* end of j-loop */

    rotate_aux_arrays (helper_arrays, length);
  } /* end of i-loop */
  /* calculate energies of 5' fragments */
  int x;
  for (x = 0; x <= length; x++) {
    count_matrix_pt.n_ij_A_e[x] = create_hashtable_list ();
  }

  int cnt1;
  for (cnt1 = 1; cnt1 <= TURN + 1; cnt1++) {
    int c_energy = 0;
    unsigned long c_count = 1;
    result_table = &count_matrix_pt.n_ij_A_e[cnt1];
    hashtable_list_add_count(result_table, c_energy, c_count);
  }

  for (j = TURN + 2; j <= length; j++) {

    /* j-1 is unpaired ... */
    matrices->f5[j] = MIN2(matrices->f5[j], matrices->f5[j - 1]);

        source_table = &count_matrix_pt.n_ij_A_e[j - 1];
        result_table = &count_matrix_pt.n_ij_A_e[j];
        for (int ei = 0; ei < source_table->length; ei++) {
          int c_energy = lookup_energy_from_index (source_table, ei);
          unsigned long c_count = lookup_count_from_index (source_table, ei);
          hashtable_list_add_count(result_table, c_energy, c_count);
        }

    /* j pairs with 1 */
    ij = fc->jindx[j] + 1;
    int type = fc->ptype[ij];
    int additional_en = 0;

    if (type) {
      if (dangles == 2)
        additional_en += E_ExtLoop (type, -1, j < length ? S1[j + 1] : -1, P);
      else
        additional_en += E_ExtLoop (type, -1, -1, P);
    }

    if (matrices->c[ij] != INF) {
      matrices->f5[j] = MIN2(matrices->f5[j], matrices->c[ij] + additional_en);

        source_table = &count_matrix_pt.n_ij_e[ij];
        result_table = &count_matrix_pt.n_ij_A_e[j];
        for (int ei = 0; ei < source_table->length; ei++) {
          int c_energy = lookup_energy_from_index (source_table, ei);
          unsigned long c_count = lookup_count_from_index (source_table, ei);
          int sum_energy = c_energy + additional_en;
          if ((sum_energy >= min_energy) && (sum_energy <= max_energy)){
            hashtable_list_add_count(result_table, sum_energy, c_count);
          }
        }
    }

    /* j pairs with some other nucleotide -> see below */
    for (i = j - TURN - 1; i > 1; i--) {
      ij = fc->jindx[j] + i;
      type = fc->ptype[ij];
      if (type) {
        if (dangles == 2)
          additional_en = E_ExtLoop (type, S1[i - 1], j < length ? S1[j + 1] : -1, P);
        else
          additional_en = E_ExtLoop (type, -1, -1, P);

        //  if (!fc->matrices->c[ij])
        //    continue;

        if (fc->matrices->f5[i - 1] != INF && fc->matrices->c[ij] != INF) {
          matrices->f5[j] = MIN2(matrices->f5[j], matrices->f5[i - 1] + matrices->c[ij] + additional_en);

        source_table = &count_matrix_pt.n_ij_e[ij];
        source_table_2 = &count_matrix_pt.n_ij_A_e[i - 1];
        result_table = &count_matrix_pt.n_ij_A_e[j];
        for (int ei_1 = 0; ei_1 < source_table->length; ei_1++) {
          int c_energy = lookup_energy_from_index (source_table, ei_1);
          unsigned long c_count = lookup_count_from_index (source_table, ei_1);
          for (int ei_2 = 0; ei_2 < source_table_2->length; ei_2++) {
              int f5_energy = lookup_energy_from_index (source_table_2, ei_2);
              unsigned long f5_count = lookup_count_from_index (source_table_2, ei_2);
              int sum_energy = c_energy + f5_energy + additional_en;
              if ((sum_energy >= min_energy) && (sum_energy <= max_energy)){
                unsigned long product_count = c_count * f5_count;
                hashtable_list_add_count(result_table, sum_energy, product_count);
              }
          }
        }
        }
      }
    }
  }

  printf ("\nF_Matrix:");
  print_array_energy_counts (count_matrix_pt.n_ij_A_e, length, min_energy, step_energy, energy_length);

  /*
   printf("\n");
   int fidx = 0;
   for(; fidx < length+1; fidx++){
   printf(" %d", matrices->f5[fidx]);
   }
   printf("\n");



   // print_matrix(f5, fc->jindx, 1, length+1);
   for(energy_index = 0; energy_index < energy_length; energy_index++){
   if(count_matrix_pt.n_ij_A_e[length]){
   int counts = count_matrix_pt.n_ij_A_e[length][energy_index];
   int energy = min_energy + energy_index * step_energy;
   printf("energy: %d, counts: %d \n", energy, counts);
   }
   }
   */

  /* clean up memory */
  free_aux_arrays (helper_arrays);

  for (i = length - turn - 1; i >= 1; i--) {
    for (j = i + turn + 1; j <= length; j++) {
      ij = indx[j] + i;

      if(count_matrix_pt.n_ij_e[ij].allocated_size > 0)
        free_hashtable_list(&count_matrix_pt.n_ij_e[ij]);
      if(count_matrix_pt.n_ij_M_e[ij].allocated_size > 0)
        free_hashtable_list(&count_matrix_pt.n_ij_M_e[ij]);
      if(count_matrix_pt.n_ij_M1_e[ij].allocated_size > 0)
        free_hashtable_list(&count_matrix_pt.n_ij_M1_e[ij]);
    }
  }

  for (i = 0; i <= length; i++) {
    free_hashtable_list(&count_matrix_pt.n_ij_A_e[i]);
  }

  free (count_matrix_pt.n_ij_e);
  free (count_matrix_pt.n_ij_A_e);
  free (count_matrix_pt.n_ij_M_e);
  free (count_matrix_pt.n_ij_M1_e);

  return 0;
}

int
main()
{
  char *sequence = vrna_read_line (stdin);
  printf ("%s\n", sequence);
  //char *sequence = "GCAACCCUUAACCCUUGGGCAAC"; //"ACGUACGUUGCAACGUACGUUGCA"; // !!! "GCACUUAACCCUUGGGCAAC";
  vrna_fold_compound_t *fc;
  vrna_md_t md;
  set_model_details (&md);
  md.uniq_ML = 1;
  md.noLP = 0;
  md.circ = 0;
//  md.dangles = 0;
//md.temperature  = 37;
  vrna_param_t *params = vrna_params (&md);
//params->model_details.circ = is_circular;
  fc = vrna_fold_compound (sequence, &(params->model_details), VRNA_OPTION_DEFAULT | VRNA_OPTION_MFE);

  if (!vrna_fold_compound_prepare (fc, VRNA_OPTION_MFE)) {
    vrna_message_warning ("vrna_mfe@mfe.c: Failed to prepare vrna_fold_compound");
  }

  int result = fill_arrays (fc);

  printf ("result: %d \n", result);

  free (sequence);
  free (params);
  vrna_fold_compound_free (fc);

}
