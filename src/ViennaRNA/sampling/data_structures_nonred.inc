#ifndef DATA_STRUCTURES_NONRED_H_
#define DATA_STRUCTURES_NONRED_H_

#ifdef VRNA_NR_SAMPLING_MPFR
#include <mpfr.h>
#endif

#define DEBUG   0

/* General mpfr funtions */

#ifdef VRNA_NR_SAMPLING_MPFR
PRIVATE int precision();      /* returns precision */


PRIVATE mpfr_rnd_t default_rnd();   /* returns default rounding mode (to nearest) */


#endif

/* defines types for different sequence stetches/loops */
typedef enum _type Type;
enum stetch_type {
  NRT_NONE_TYPE       = 0,    /* nonetype: both 0 (root) */
  NRT_HAIRPIN         = 1,    /* hairpin: both 0 (unused) */
  NRT_IT_LOOP         = 2,    /* internal loop: (i', j') - closed pair */
  NRT_MB_LOOP         = 3,    /* multibranch loop: both 0 (unused) */
  NRT_EXT_LOOP        = 4,    /* external loop: (i, j) */
  NRT_UNPAIRED_SG     = 5,    /* unpaired positions: (first unpaired position - end known of stretch) */
  NRT_QM_UNPAIR       = 6,    /* QM unpaired base at 3' end : (j - 1, j) */
  NRT_QM_BRANCH       = 7,    /* QM split on K with pairing within [i..k-1] : (k, 0)  */
  NRT_QM_NOBRANCH     = 8,    /* QM split on K with no pairing within [i..k-1] : (k, 0) */
  NRT_QM2_BRANCH      = 9,    /* QM2 split on K with at least one more branch within [i..k-1] : (k, 0)  */
  NRT_QM2_UNPAIR      = 10,   /* QM2 split with unpaired bases on 3' end, i.e. j (i, j) */
  NRT_QM1_NEW_BRANCH  = 11,   /* branch in multiloop in qm1_new: (i,j) */
  NRT_EXT_HAIRPIN     = 12,   /* circRNA: ext-loop hairpin: both 0 (unused) */
  NRT_EXT_IT_LOOP     = 13,   /* circRNA: internal loop: (i', j') - closed pair */
  NRT_EXT_MB_LOOP     = 14,   /* circRNA: multibranch loop: both 0 (unused) */
};

#ifdef VRNA_NR_SAMPLING_HASH
/******************************************/
/*        version with hash table         */
/******************************************/

/* explanation: each node (tr_node) holds a pointer to a start of a linked
 * list to the next node, link to previous node, Boltzmann's factor
 * for given flat structure and itentifier, which corresponds to an id of
 * chosen flat structure
 */

/* define node of a tree */
typedef struct tr_node tr_node;


/* define elements of hashtable */
typedef struct hash_node hash_node; /* list of arrays sorting triplets for every hashvalue*/

struct tr_node {
  int       type;
  int       loop_spec_1;
  int       loop_spec_2;
  int       seqlen;
  tr_node   *parent;
  tr_node   *child;           /* used when this node has at most one child */
  hash_node *hash_tab;        /* used when there is at least two children */
  tr_node   *next_in_hash;    /* used for node that is next in hash table */
  int       hash_value;       /* hash value computed for this node */
#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_t    weight;
  mpfr_t    max_weight;   /* maximum allowed weight (maximum of partition function) */
#else
  double    weight;
  double    max_weight;       /* maximum allowed weight (maximum of partition function) */
#endif
  int       created_recently; /* 1 if was created during last iteration, otherwise 0 */
};

struct hash_node {
  tr_node **array;
  int     hash_elemts;
  int     hash_size;
};

/* tree + hash functions */
/** @brief creates a root of datastructure tree (hash table version) **/
PRIVATE tr_node *create_root(int    seqlen,
                             double max_weight);


/** @brief returns a weigh of node (type, loop_spec_1, loop_spec_2) if child of last_node, otherwise returns 0.0 **/
PRIVATE double tr_node_weight(tr_node *last_node,
                              int     type,
                              int     loop_spec_1,
                              int     loop_spec_2);


/** @brief sums weight of all children of par_node and returns it **/
PRIVATE double total_weight_par(tr_node *par_node);


/** @brief sums weight of all children of par_node with certain type and returns it **/
PRIVATE double total_weight_par_type(int      type,
                                     tr_node  *par_node);


/** @brief creates node (type, loop_spec_1, loop_spec_2) if not existing and returns pointer to it,
 * or returns pointer to exisiting case **/
PRIVATE tr_node *add_if_nexists(int     type,
                                int     loop_spec_1,
                                int     loop_spec_2,
                                tr_node *par_node,
                                double  max_weight);


/** @brief traces back from leaf to root while updating weights of leaf to all nodes in path,
 *  returns pointer to root **/
PRIVATE tr_node *traceback_to_root(tr_node  *leaf,
                                   double   struct_weight,
                                   int      *is_dup,
                                   int      *pf_overflow);


/** @brief destructor **/
PRIVATE void free_all_nr(tr_node *root);


#else
/******************************************/
/*       version with linked lists        */
/******************************************/

typedef struct tllr_node tllr_node;

struct tllr_node {
  int       type;
  int       loop_spec_1;
  int       loop_spec_2;
  tllr_node *parent;          /* vertical chaining - ancestor */
  tllr_node *head;            /* vertical chaining - successor */
  tllr_node *next_node;       /*horizontal chaining - linked list */
#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_t    weight;
  mpfr_t    max_weight;   /* maximum allowed weight (maximum of partition function) */
#else
  double    weight;
  double    max_weight;       /* maximum allowed weight (maximum of partition function) */
#endif
  int       created_recently; /* 1 if was created during last iteration, otherwise 0 */
};


/* memory object for non-redundant sampling approach using linked lists*/
typedef struct nr_memory nr_memory;

struct nr_memory {
  tllr_node *nr_memory_allocated;
  int       memory_index;
  size_t    tllr_node_size;
  size_t    block_size;   /* block size */
  nr_memory *prev_block;  /* stores previous block */
};

/* creates an object nr_memory that pre-allocates a block of memory for memory allocation */
PRIVATE nr_memory *create_nr_memory(size_t    node_size,
                                    size_t    block_size,
                                    nr_memory *prev_memory);


/* tree + linked list functions */
/** @brief creates a root of datastructure tree (linked list version) **/
PRIVATE tllr_node *create_ll_root(struct nr_memory  **memory_dat,
                                  double            max_weight);


/** resets cursor to current_node and start of linked list **/
PRIVATE void reset_cursor(tllr_node **memorized_node_prev,
                          tllr_node **memorized_node_cur,
                          tllr_node *current_node);


/** @brief moves cursor to next node if current_node is identical to one in loop, otherwise does nothing **/
PRIVATE void advance_cursor(tllr_node **memorized_node_prev,
                            tllr_node **memorized_node_cur,
                            int       type,
                            int       loop_spec_1,
                            int       loop_spec_2);


/** @brief returns a weigh of node (type, loop_spec_1, loop_spec_2) if child of last_node, otherwise returns 0.0 **/
PRIVATE double get_weight(tllr_node *memorized_node_cur,
                          int       type,
                          int       loop_spec_1,
                          int       loop_spec_2);


/** @brief sums weight of all children of par_node and returns it **/
PRIVATE double get_weight_all(tllr_node *par_node);


/** @brief sums weight of all children of par_node with certain type and returns it **/
PRIVATE double get_weight_type_spec(int       type,
                                    tllr_node *par_node);


/** @brief creates node (type, loop_spec_1, loop_spec_2) if not existing and returns pointer to it,
 * or returns pointer to exisiting case **/
PRIVATE tllr_node *add_if_nexists_ll(struct nr_memory **memory_dat,
                                     int              type,
                                     int              loop_spec_1,
                                     int              loop_spec_2,
                                     tllr_node        *memorized_node_prev,
                                     tllr_node        *memorized_node_cur,
                                     tllr_node        *parent_node,
                                     double           max_weight);


/** @brief traces back from leaf to root while updating weights of leaf to all nodes in path,
 *  returns pointer to root **/
PRIVATE tllr_node *traceback_to_ll_root(tllr_node *leaf,
                                        double    weight,
                                        int       *is_dup,
                                        int       *pf_overflow);


/** @brief destructor **/
PRIVATE void free_all_nrll(nr_memory **memory_dat);


#endif


/**************************************************/
/*            mpfr related functions              */
/**************************************************/

#ifdef VRNA_NR_SAMPLING_MPFR
PRIVATE int
precision()
{
  /* returns precision */
  return 128;
}


PRIVATE mpfr_rnd_t
default_rnd()
{
  /* returns default rounding mode (to nearest) */
  return mpfr_get_default_rounding_mode();
}


#endif

/**************************************************/
/*                node weight getter              */
/**************************************************/

#ifdef VRNA_NR_SAMPLING_HASH
PRIVATE double
return_node_weight(tr_node *node)
{
#ifdef VRNA_NR_SAMPLING_MPFR
  return mpfr_get_d(node->weight, default_rnd());
#else
  return node->weight;
#endif
}


#else
PRIVATE double
return_node_weight(tllr_node *node)
{
#ifdef VRNA_NR_SAMPLING_MPFR
  return mpfr_get_d(node->weight, default_rnd());
#else
  return node->weight;
#endif
}


#endif

/**************************************************/
/*       data structure (tree + hash table)       */
/**************************************************/
/* Contains
 *
 *
 *
 * Creates a structure memorizing the weights of already used solutions for the
 * non-redundant sampling. It uses tree structure where each node contains a hash
 * table.
 */


#ifdef VRNA_NR_SAMPLING_HASH


PRIVATE int
pow_of_2_modulus(int  numerator,
                 int  denominator)
{
  return numerator & (denominator - 1);
}


/* creates new tree node with undefined weight with parent node */
PRIVATE tr_node *
create_tr_node(int      type,
               int      loop_spec_1,
               int      loop_spec_2,
               int      seqlen,
               tr_node  *parent,
               double   max_weight)
{
  tr_node *new_tr_node = (tr_node *)vrna_alloc(sizeof(tr_node));

  new_tr_node->type = type;
  /* Types and properties specific to loops:
   * type 0 : nonetype: both 0 (root)
   * type 1 : hairpin: both 0 (unused)
   * type 2 : internal loop: (i', j') - closed pair
   * type 3 : multiloop: K (second is 0 - unused)
   * type 4 : external loop: (i, j)
   * type 5 : unpaired positions: (first unpaired position - end known of stretch)
   * type 6 : special type qm1: (i, j) (used in multiloops)
   * type 7 : branch in multiloop in qm1: (i,j)
   */
  new_tr_node->loop_spec_1  = loop_spec_1;
  new_tr_node->loop_spec_2  = loop_spec_2;
  new_tr_node->seqlen       = seqlen;
  new_tr_node->parent       = parent;
  new_tr_node->child        = NULL;
  new_tr_node->hash_tab     = NULL;
  new_tr_node->next_in_hash = NULL;
  new_tr_node->hash_value   = 0;
#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_init2(new_tr_node->weight, precision());
  mpfr_set_d(new_tr_node->weight, 0., default_rnd());
  mpfr_init2(new_tr_node->max_weight, precision());
  mpfr_set_d(new_tr_node->max_weight, max_weight, default_rnd());
#else
  new_tr_node->weight     = 0;
  new_tr_node->max_weight = max_weight;
#endif
  new_tr_node->created_recently = 1;
  return new_tr_node;
}


/* compares hash children values with actual value in parent */
PRIVATE void
compare_parent_children_weight(tr_node *parent)
{
  tr_node *node_t;
  int     i;

#ifdef  VRNA_NR_SAMPLING_MPFR
  mpfr_t  total;
  mpfr_init2(total, precision());
  mpfr_set_d(total, 0., default_rnd());
#else
  double  total = 0.;
#endif

  if (!parent->hash_tab) {
    if (parent->child) {
#ifdef VRNA_NR_SAMPLING_MPFR
      mpfr_add(total, total, parent->child->weight, default_rnd());
#else
      total += parent->child->weight;
#endif
    }
  } else {
    for (i = 0; i < parent->hash_tab->hash_size; i++) {
      node_t = parent->hash_tab->array[i];
      while (node_t) {
#ifdef VRNA_NR_SAMPLING_MPFR
        mpfr_add(total, total, node_t->weight, default_rnd());
#else
        total += node_t->weight;
#endif
        node_t = node_t->next_in_hash;
      }
    }
  }

#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_clear(total);
#endif
}


/* creates hash entry point (its start) */
PRIVATE hash_node *
init_hash()
{
  int       i;
  hash_node *new_hash = (hash_node *)vrna_alloc(sizeof(hash_node));

  new_hash->hash_elemts = 0;
  new_hash->hash_size   = 4; /* will be expanded when necessary */
  new_hash->array       = (tr_node **)vrna_alloc(sizeof(tr_node *) * (new_hash->hash_size));
  for (i = 0; i < new_hash->hash_size; i++)
    new_hash->array[i] = NULL;
  return new_hash;
}


/* creates root (start of a tree) */
PRIVATE tr_node *
create_root(int     seqlen,
            double  max_weight)
{
  tr_node *root = create_tr_node(NRT_NONE_TYPE, 0, 0, seqlen, NULL, max_weight);

  return root;
}


/* returns hash key for given triple of values */
PRIVATE int
hash_fci(int  type,
         int  loop_spec_1,
         int  loop_spec_2,
         int  seqlen)
{
  int hash_f = 3 * type + 11 * loop_spec_1 + 50021 * loop_spec_2;

  /* int hash_f = type + 9*loop_spec_1 + 9*(seqlen+1)*loop_spec_2; */
  return hash_f;
}


/* Expands hash when its load becomes higher than 50% number and relinks elements.
 * New size is lowest prime number higher than twice the size.
 * */
PRIVATE void
extend_hash(hash_node *hash_to_res)
{
  /* extends table */
  tr_node *next;
  tr_node *tmp;
  int     hash_val, hash_val_mod;
  int     hash_size_tmp = hash_to_res->hash_size;

  hash_to_res->hash_size  *= 2;
  hash_to_res->array      =
    (tr_node **)realloc(hash_to_res->array, sizeof(tr_node *) * (hash_to_res->hash_size));
  for (int i = hash_size_tmp; i < hash_to_res->hash_size; i++)
    hash_to_res->array[i] = NULL;
  /* relinks all elements	within table
   * this table points to beginning of each chain within hash table */
  tr_node **hash_tmp = (tr_node **)vrna_alloc(sizeof(tr_node *) * (hash_to_res->hash_size));
  for (int i = 0; i < hash_to_res->hash_size; i++)
    hash_tmp[i] = NULL;

  /* this table points at the current end of eaxh hash table */
  tr_node **hash_cur_end = (tr_node **)vrna_alloc(sizeof(tr_node *) * (hash_to_res->hash_size));
  for (int i = 0; i < hash_to_res->hash_size; i++)
    hash_cur_end[i] = NULL;
  for (int i = 0; i < hash_size_tmp; i++) {
    next = hash_to_res->array[i];
    while (next) {
      hash_val      = hash_fci(next->type, next->loop_spec_1, next->loop_spec_2, next->seqlen);
      hash_val_mod  = pow_of_2_modulus(hash_val, hash_to_res->hash_size);
      if (!hash_tmp[hash_val_mod]) {
        hash_tmp[hash_val_mod]      = next;
        hash_cur_end[hash_val_mod]  = next;
      } else {
        hash_cur_end[hash_val_mod]->next_in_hash  = next;
        hash_cur_end[hash_val_mod]                = next;
      }

      tmp               = next;
      next              = next->next_in_hash;
      tmp->next_in_hash = NULL;
    }
  }
  for (int i = 0; i < hash_to_res->hash_size; i++)
    hash_to_res->array[i] = hash_tmp[i];
  free(hash_tmp);
  free(hash_cur_end);
}


/* inserting new value into corresponding hash */
PRIVATE void
insert_val(hash_node  *hash_t,
           int        type,
           int        loop_spec_1,
           int        loop_spec_2,
           tr_node    *inserted_node)
{
  int     hash_val_mod;

  tr_node *ptr;

  /* where to insert */
  inserted_node->hash_value = hash_fci(type, loop_spec_1, loop_spec_2, inserted_node->seqlen);
  hash_val_mod              = pow_of_2_modulus(inserted_node->hash_value, hash_t->hash_size);


  if (!hash_t->array[hash_val_mod]) {
    /* head case */
    hash_t->array[hash_val_mod] = inserted_node;
    return;
  }

  ptr = hash_t->array[hash_val_mod];
  while (ptr->next_in_hash) /* inside list */
    ptr = ptr->next_in_hash;
  ptr->next_in_hash = inserted_node;
}


/* creates and inserts a new tr_node (along with creating hash) into tree */
PRIVATE tr_node *
insert_tr_node(int      type,
               int      loop_spec_1,
               int      loop_spec_2,
               tr_node  *last_node,
               double   max_weight)
{
  /* forge new tr_node */
  tr_node *new_tr_node =
    create_tr_node(type, loop_spec_1, loop_spec_2, last_node->seqlen, last_node, max_weight);

  if (!last_node->child) {
    /* empty so no child yet - insert it directly */
    last_node->child = new_tr_node;
  } else {
    if (!last_node->hash_tab) {
      /* hash not defined yet (we have exactly two values) */
      last_node->hash_tab = init_hash(last_node->seqlen);
      insert_val(last_node->hash_tab,
                 last_node->child->type,
                 last_node->child->loop_spec_1,
                 last_node->child->loop_spec_2,
                 last_node->child);
      last_node->hash_tab->hash_elemts++;
    }

    /* insert it into hashtable of previous node */
    insert_val(last_node->hash_tab, type, loop_spec_1, loop_spec_2, new_tr_node);
    last_node->hash_tab->hash_elemts++;
    if (((double)last_node->hash_tab->hash_elemts / last_node->hash_tab->hash_size) > 0.5)
      extend_hash(last_node->hash_tab);
  }

  return new_tr_node;
}


/* tracks whether there is node in hash with given values */
PRIVATE tr_node *
access_val(tr_node  *parent,
           int      type,
           int      loop_spec_1,
           int      loop_spec_2)
{
  tr_node *ptr = NULL;

  if (!parent->hash_tab) {
    /* hash does not exist (at most one child) */
    if (parent->child) {
      /* child exists (exactly one child) */
      if ((parent->child->type == type) && (parent->child->loop_spec_1 == loop_spec_1)
          && (parent->child->loop_spec_2 == loop_spec_2))
        return parent->child;
    }
  } else {
    /* hash exists -> 2+ children */
    int hash_f      = hash_fci(type, loop_spec_1, loop_spec_2, parent->seqlen);
    int hash_f_mod  = pow_of_2_modulus(hash_f, parent->hash_tab->hash_size);
    ptr = parent->hash_tab->array[hash_f_mod];
    while (ptr) {
      if (ptr->hash_value == hash_f)
        break;

      ptr = ptr->next_in_hash;
    }
  }

  if (ptr)
    return ptr;

  return NULL;
}


/* creates node with properties that don't already exist for given node */
PRIVATE tr_node *
add_if_nexists(int      type,
               int      loop_spec_1,
               int      loop_spec_2,
               tr_node  *last_node,
               double   max_weight)
{
  tr_node *new_tr_node = access_val(last_node, type, loop_spec_1, loop_spec_2);

  if (!new_tr_node)
    new_tr_node = insert_tr_node(type, loop_spec_1, loop_spec_2, last_node, max_weight);

  return new_tr_node;
}


/* returns weight of node if it exists; otherwise returns 0 */
PRIVATE double
tr_node_weight(tr_node  *last_node,
               int      type,
               int      loop_spec_1,
               int      loop_spec_2)
{
  tr_node *next_node = access_val(last_node, type, loop_spec_1, loop_spec_2);

  if (next_node) {
#ifdef VRNA_NR_SAMPLING_MPFR
    return mpfr_get_d(next_node->weight, default_rnd());
#else
    return next_node->weight;
#endif
  }

  return 0.;
}


/* returns sum of all weights for given parent */
PRIVATE double
total_weight_par(tr_node *last_node)
{
  if ((!last_node->hash_tab) && (!last_node->child))
    return 0.;

#ifdef VRNA_NR_SAMPLING_MPFR
  return mpfr_get_d(last_node->weight, default_rnd());
#else
  return last_node->weight;
#endif
}


/* returns sum of all weights for given parent and for given type*/
PRIVATE double
total_weight_par_type(int     type,
                      tr_node *last_node)
{
  tr_node *ptr = NULL;
  int     i;

#ifdef VRNA_NR_SAMPLING_MPFR
  double  weight_sum_d;
  mpfr_t  weight_sum;
  mpfr_init2(weight_sum, precision());
  mpfr_set_d(weight_sum, 0., default_rnd());
#else
  double  weight_sum = 0.;
#endif


  if (!last_node->hash_tab) {
    if (last_node->child)
      if (last_node->child->type == type) {
#ifdef VRNA_NR_SAMPLING_MPFR
        mpfr_clear(weight_sum);
        return mpfr_get_d(last_node->child->weight, default_rnd());
#else
        return last_node->child->weight;
#endif
      }
  } else {
    for (i = 0; i < last_node->hash_tab->hash_size; i++) {
      ptr = last_node->hash_tab->array[i];
      while (ptr) {
        if (ptr->type == type) {
#ifdef VRNA_NR_SAMPLING_MPFR
          mpfr_add(weight_sum, weight_sum, ptr->weight, default_rnd());
#else
          weight_sum += ptr->weight;
#endif
        }

        ptr = ptr->next_in_hash;
      }
    }
  }

#ifdef VRNA_NR_SAMPLING_MPFR
  weight_sum_d = mpfr_get_d(weight_sum, default_rnd());
  mpfr_clear(weight_sum);
  return weight_sum_d;
#else
  return weight_sum;
#endif
}


/* updates weight of a given node */
PRIVATE int
update_weight(tr_node *node,
              double  weight)
{
#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_t intermediate;
  mpfr_init2(intermediate, precision());
  mpfr_add_d(intermediate, node->weight, weight, default_rnd());
  mpfr_sub(intermediate, node->max_weight, intermediate, default_rnd());

  if (mpfr_cmp_d(intermediate, -1E-14) < 0) {
    mpfr_clear(intermediate);
    return 1;
  } else {
    mpfr_clear(intermediate);
    mpfr_add_d(node->weight, node->weight, weight, default_rnd());
  }

#else
  if (node->max_weight - (node->weight + weight) < -(1E-14))
    return 1;
  else
    node->weight += weight;

#endif
  return 0;
}


/* tracebacks to root while updating values for each node passed through
 * - also verifies unicity (at least one node differs) */
PRIVATE tr_node *
traceback_to_root(tr_node *leaf,
                  double  struct_weight,
                  int     *is_dup,
                  int     *pf_overflow)
{
  *pf_overflow = update_weight(leaf, struct_weight);
  if (leaf->created_recently) {
    leaf->created_recently  = 0;
    *is_dup                 = 0;
  }

  while (leaf->parent) {
    *pf_overflow = update_weight(leaf->parent, struct_weight);
    if (leaf->parent->created_recently) {
      leaf->parent->created_recently  = 0;
      *is_dup                         = 0;
    }

    leaf = leaf->parent;
  }
  return leaf;
}


/* destructor */
PRIVATE void
free_all_nr(tr_node *root)
{
  if (!root->hash_tab) {
    if (root->child)
      free_all_nr(root->child);
  } else {
    for (int i = 0; i < root->hash_tab->hash_size; i++) {
      tr_node *ptr = root->hash_tab->array[i];
      tr_node *tmp;
      while (ptr) {
        tmp = ptr->next_in_hash;
        free_all_nr(ptr);
        ptr = tmp;
      }
    }
    free(root->hash_tab->array);
    free(root->hash_tab);
  }

#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_clear(root->weight);
  mpfr_clear(root->max_weight);
#endif
  free(root);
}


#endif


/*********************************************************/
/*       data structure (linked_list + hash table)       */
/*********************************************************/

#ifndef VRNA_NR_SAMPLING_HASH
/* allocates a  nr_memory object -  pre-allocator for tllr_nodes */
PRIVATE nr_memory *
create_nr_memory(size_t     tllr_node_size,
                 size_t     block_size,
                 nr_memory  *prev_block)
{
  struct nr_memory *memory_dat = vrna_alloc(sizeof(nr_memory));

  memory_dat->nr_memory_allocated = vrna_alloc(block_size);
  memory_dat->memory_index        = 0;
  memory_dat->tllr_node_size      = tllr_node_size;
  memory_dat->block_size          = block_size;
  memory_dat->prev_block          = prev_block;

  return memory_dat;
}


/* This creates structure that uses linked list instead of hash. The thought behind this is
 * the order of investigated nodes is always the same so we can add them to specific place.
 * It is thus a bit faster.
 */
PRIVATE tllr_node *
create_tllr_node(struct nr_memory **memory_dat,
                 int              type,
                 int              loop_spec_1,
                 int              loop_spec_2,
                 tllr_node        *parent,
                 double           max_weight)
{
  tllr_node *new_tllr_node;

  if ((*memory_dat)->tllr_node_size * ((*memory_dat)->memory_index + 1) >
      (*memory_dat)->block_size) {
    struct nr_memory *memory_dat_tmp = create_nr_memory((*memory_dat)->tllr_node_size,
                                                        (*memory_dat)->block_size,
                                                        *memory_dat);
    *memory_dat   = memory_dat_tmp;
    new_tllr_node = &((*memory_dat)->nr_memory_allocated[(*memory_dat)->memory_index]);
  } else {
    new_tllr_node = &((*memory_dat)->nr_memory_allocated[(*memory_dat)->memory_index]);
  }

  new_tllr_node->type = type;
  /* Types and properties specific to loops:
   * type 0 : nonetype: both 0 (root)
   * type 1 : hairpin: both 0 (unused)
   * type 2 : internal loop: (i', j') - closed pair
   * type 3 : multiloop: K (second is 0 - unused)
   * type 4 : external loop: (i, j)
   * type 5 : unpaired positions: (first unpaired position - end known of stretch)
   * type 6 : special type qm1: (i, j) (used in multiloops)
   * type 7 : branch in multiloop in qm1: (i,j)
   */
  new_tllr_node->loop_spec_1  = loop_spec_1;
  new_tllr_node->loop_spec_2  = loop_spec_2;
  new_tllr_node->parent       = parent;
  new_tllr_node->next_node    = NULL;
  new_tllr_node->head         = NULL;
#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_init2(new_tllr_node->weight, precision());
  mpfr_set_d(new_tllr_node->weight, 0., default_rnd());
  mpfr_init2(new_tllr_node->max_weight, precision());
  mpfr_set_d(new_tllr_node->max_weight, max_weight, default_rnd());
#else
  new_tllr_node->weight     = 0;
  new_tllr_node->max_weight = max_weight;
#endif
  new_tllr_node->created_recently = 1;

  (*memory_dat)->memory_index++;
  return new_tllr_node;
}


/* compares hash children values with actual value in parent */
#if DEBUG
PRIVATE void
compare_parent_children_weight_tr(tllr_node *parent)
{
  tllr_node *node_t;

#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_t    total;
  mpfr_init2(total, precision());
  mpfr_set_d(total, 0., default_rnd());
#else
  double    total = 0.;
#endif

  node_t = parent->head;
  while (node_t) {
#ifdef VRNA_NR_SAMPLING_MPFR
    mpfr_add(total, total, node_t->weight, default_rnd());
#else
    total += node_t->weight;
#endif
    node_t = node_t->next_node;
  }
#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_clear(total);
#endif
}


#endif

/* creates root (start of a tree) */
PRIVATE tllr_node *
create_ll_root(struct nr_memory **memory_dat,
               double           max_weight)
{
  tllr_node *root = create_tllr_node(memory_dat, NRT_NONE_TYPE, 0, 0, NULL, max_weight);

  return root;
}


/* inserts a tllr_node before 'next_node' and after previous node
 * Node cursor hols previous and current ll_node */
PRIVATE tllr_node *
insert_tllr_node(struct nr_memory **memory_dat,
                 tllr_node        *memorized_node_prev,
                 tllr_node        *memorized_node_cur,
                 int              type,
                 int              loop_spec_1,
                 int              loop_spec_2,
                 tllr_node        *parent_node,
                 double           max_weight)
{
  tllr_node *new_node = create_tllr_node(memory_dat,
                                         type,
                                         loop_spec_1,
                                         loop_spec_2,
                                         parent_node,
                                         max_weight);

  if (!memorized_node_prev) /* first node to be inserted */
    parent_node->head = new_node;
  else
    memorized_node_prev->next_node = new_node;

  new_node->next_node = memorized_node_cur;
  return new_node;
}


/* resets cursor to beginning of loop*/
PRIVATE void
reset_cursor(tllr_node  **memorized_node_prev,
             tllr_node  **memorized_node_cur,
             tllr_node  *current_node)
{
  (*memorized_node_prev)  = NULL;
  (*memorized_node_cur)   = current_node->head;
}


/* advances pointer in loop if the identifier coincide with current pointer and returns weight */
PRIVATE inline void
advance_cursor(tllr_node  **memorized_node_prev,
               tllr_node  **memorized_node_cur,
               int        type,
               int        loop_spec_1,
               int        loop_spec_2)
{
  if (*memorized_node_cur) {
    if ((*memorized_node_cur)->type == type
        && (*memorized_node_cur)->loop_spec_1 == loop_spec_1
        && (*memorized_node_cur)->loop_spec_2 == loop_spec_2) {
      (*memorized_node_prev)  = (*memorized_node_cur);
      (*memorized_node_cur)   = (*memorized_node_cur)->next_node;
    }
  }
}


/* gets weight of actual node */
PRIVATE inline double
get_weight(tllr_node  *memorized_node_cur,
           int        type,
           int        loop_spec_1,
           int        loop_spec_2)
{
  double weight = 0;

  if (memorized_node_cur) {
    if (memorized_node_cur->type == type
        && memorized_node_cur->loop_spec_1 == loop_spec_1
        && memorized_node_cur->loop_spec_2 == loop_spec_2)
#ifdef VRNA_NR_SAMPLING_MPFR
      weight = mpfr_get_d(memorized_node_cur->weight, default_rnd());

#else
      weight = memorized_node_cur->weight;
#endif
  }

  return weight;
}


/* get weight of all child nodes */
PRIVATE double
get_weight_all(tllr_node *last_node)
{
  if (!last_node->head)
    return 0;

#ifdef VRNA_NR_SAMPLING_MPFR
  return mpfr_get_d(last_node->weight, default_rnd());
#else
  return last_node->weight;
#endif
}


/* get weight of all child nodes of certain type */
PRIVATE double
get_weight_type_spec(int        type,
                     tllr_node  *last_node)
{
  /* double    weight_total  = 0; */

#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_t    weight_total;
  double    weight_total_d;
  mpfr_init2(weight_total, precision());
  mpfr_set_d(weight_total, 0., default_rnd());
#else
  double    weight_total = 0.;
#endif

  tllr_node *ptr = last_node->head;

  while (ptr) {
    if (ptr->type == type) {
#ifdef VRNA_NR_SAMPLING_MPFR
      mpfr_add(weight_total, weight_total, ptr->weight, default_rnd());
#else
      weight_total += ptr->weight;
#endif
    }

    ptr = ptr->next_node;
  }

#ifdef VRNA_NR_SAMPLING_MPFR
  weight_total_d = mpfr_get_d(weight_total, default_rnd());
  mpfr_clear(weight_total);
  return weight_total_d;
#else
  return weight_total;
#endif
}


/* adds node if the current one isn't the one we want */
PRIVATE inline tllr_node *
add_if_nexists_ll(struct nr_memory  **memory_dat,
                  int               type,
                  int               loop_spec_1,
                  int               loop_spec_2,
                  tllr_node         *memorized_node_prev,
                  tllr_node         *memorized_node_cur,
                  tllr_node         *parent_node,
                  double            max_weight)
{
  tllr_node *returned_node;

  if (memorized_node_cur) {
    if (memorized_node_cur->type == type
        && memorized_node_cur->loop_spec_1 == loop_spec_1
        && memorized_node_cur->loop_spec_2 == loop_spec_2)
      returned_node = memorized_node_cur;
    else
      returned_node = insert_tllr_node(memory_dat,
                                       memorized_node_prev,
                                       memorized_node_cur,
                                       type,
                                       loop_spec_1,
                                       loop_spec_2,
                                       parent_node,
                                       max_weight);
  } else {
    returned_node = insert_tllr_node(memory_dat,
                                     memorized_node_prev,
                                     memorized_node_cur,
                                     type,
                                     loop_spec_1,
                                     loop_spec_2,
                                     parent_node,
                                     max_weight);
  }

  return returned_node;
}


/* updates weight of a given node */
PRIVATE int
update_weight_ll(tllr_node  *node,
                 double     weight)
{
#ifdef VRNA_NR_SAMPLING_MPFR
  mpfr_t intermediate;
  mpfr_init2(intermediate, precision());
  mpfr_add_d(intermediate, node->weight, weight, default_rnd());
  mpfr_sub(intermediate, node->max_weight, intermediate, default_rnd());

  /* if(node->max_weight - (node->weight + weight) < -(1E-14)){ */
  if (mpfr_cmp_d(intermediate, -1E-14) < 0) {
    mpfr_clear(intermediate);
    return 1;
  } else {
    mpfr_clear(intermediate);
    mpfr_add_d(node->weight, node->weight, weight, default_rnd());
  }

#else
  if (node->max_weight - (node->weight + weight) < -(1E-14))
    return 1;
  else
    node->weight += weight;

#endif
  return 0;
}


/* tracebacks to root while updating values for each node passed through
 * - also verifies unicity (at least one node differs) */
PRIVATE tllr_node *
traceback_to_ll_root(tllr_node  *leaf,
                     double     weight,
                     int        *is_dup,
                     int        *pf_overflow)
{
  *pf_overflow = update_weight_ll(leaf, weight);
  if (leaf->created_recently) {
    /* check whether the last sequence is not a duplicate */
    leaf->created_recently  = 0;
    *is_dup                 = 0;
  }

  while (leaf->parent) {
    *pf_overflow = update_weight_ll(leaf->parent, weight);
    if (leaf->parent->created_recently) {
      leaf->parent->created_recently  = 0;
      *is_dup                         = 0;
    }

    leaf = leaf->parent;
  }
  return leaf;
}


/* destructor */
PRIVATE void
free_all_nrll(struct nr_memory **memory_dat)
{
  int i;

  if (memory_dat) {
    nr_memory *memory_block = *memory_dat;
    nr_memory *memory_block_next;
    while (memory_block) {
      memory_block_next = memory_block->prev_block;
#ifdef VRNA_NR_SAMPLING_MPFR
      for (i = 0; i < memory_block->memory_index; i++) {
        mpfr_clear(memory_block->nr_memory_allocated[i].weight);
        mpfr_clear(memory_block->nr_memory_allocated[i].max_weight);
      }
#endif
      free(memory_block->nr_memory_allocated);
      free(memory_block);
      memory_block = memory_block_next;
    }
  }
}


#endif


#endif
