/**************************************************/
/*       data structure (tree + hash table)       */
/**************************************************/
/* Contains
 * 


/* Creates a structure memorizing the weights of already used solutions for the 
 * non-redundant sampling. It uses tree structure where each node contains a hash
 * table. 
 */

#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "data_structures_nonred.h"

int pow_of_2_modulus(int numerator, int denominator){
	return (numerator & (denominator-1));
}

/* creates new tree node with undefined weight with parent node */
tr_node * create_tr_node(int type, int loop_spec_1, int loop_spec_2, int seqlen, tr_node * parent){
  tr_node * new_tr_node = (tr_node *) malloc(sizeof(tr_node));
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
  new_tr_node->loop_spec_1 = loop_spec_1;
  new_tr_node->loop_spec_2 = loop_spec_2;
  new_tr_node->seqlen = seqlen;
  new_tr_node->parent = parent;
  new_tr_node->child = NULL;
  new_tr_node->hash_tab = NULL;
  new_tr_node->next_in_hash = NULL;
  new_tr_node->hash_value = 0;
  new_tr_node->weight = 0;
  return new_tr_node;
}

/* creates hash entry point (its start) */
hash_node * init_hash(){
	int i;
	hash_node* new_hash = (hash_node *) malloc(sizeof(hash_node));
	new_hash->hash_elemts = 0;
	new_hash->hash_size = 4; /* will be expanded when necessary */
	new_hash->array = (tr_node **) malloc(sizeof(tr_node*) * (new_hash->hash_size));
	for(i=0; i<new_hash->hash_size; i++){
		new_hash->array[i]=NULL;
	}
	return new_hash;
}

/* creates root (start of a tree) */
tr_node * create_root(int seqlen){
  tr_node * root = create_tr_node(NRT_NONE_TYPE, 0, 0, seqlen, NULL);
  return root;  	
}


/* returns hash key for given triple of values */
int hash_fci(int type, int loop_spec_1, int loop_spec_2, int seqlen){
	int hash_f = 3*type + 11*loop_spec_1 + 50021*loop_spec_2;
	//int hash_f = type + 9*loop_spec_1 + 9*(seqlen+1)*loop_spec_2;
	return hash_f;
}


/* Expands hash when its load becomes higher than 50% number and relinks elements.
 * New size is lowest prime number higher than twice the size.
 * */
void extend_hash(hash_node * hash_to_res){
	/* extends table */
	tr_node * next;
	tr_node * link;
	tr_node * tmp;
	int hash_val, hash_val_mod;
	int hash_size_tmp = hash_to_res->hash_size;	
	hash_to_res->hash_size *= 2;
	hash_to_res->array = (tr_node **)realloc(\
	 hash_to_res->array, sizeof(tr_node*) * (hash_to_res->hash_size));
	for(int i=hash_size_tmp; i<hash_to_res->hash_size; i++){
		hash_to_res->array[i]=NULL;
	}
	/* relinks all elements	within table
	 * this table points to beginning of each chain within hash table */
	tr_node ** hash_tmp = (tr_node**)malloc(sizeof(tr_node*) * (hash_to_res->hash_size));
	for(int i=0; i<hash_to_res->hash_size; i++){
		hash_tmp[i] = NULL;
	}
	
	/* this table points at the current end of eaxh hash table */
	tr_node ** hash_cur_end = (tr_node**)malloc(sizeof(tr_node*) * (hash_to_res->hash_size));
	for(int i=0; i<hash_to_res->hash_size; i++){
		hash_cur_end[i] = NULL;
	}
	for(int i=0; i<hash_size_tmp; i++){
		next = hash_to_res->array[i];
		while(next){
			hash_val = hash_fci(next->type, next->loop_spec_1, next->loop_spec_2, next->seqlen);
			hash_val_mod = pow_of_2_modulus(hash_val, hash_to_res->hash_size);
			if(!hash_tmp[hash_val_mod]){
				hash_tmp[hash_val_mod] = next;
				hash_cur_end[hash_val_mod] = next;
			}else{
				hash_cur_end[hash_val_mod]->next_in_hash = next;
				hash_cur_end[hash_val_mod] = next;
			}
			tmp = next;
			next = next->next_in_hash;
			tmp->next_in_hash = NULL;
		}
	}
	for(int i=0; i<hash_to_res->hash_size; i++){
		hash_to_res->array[i] = hash_tmp[i];	
	} 
	free(hash_tmp);
	free(hash_cur_end);
}


/* inserting new value into corresponding hash */
void insert_val(hash_node * hash_t, int type, int loop_spec_1, int loop_spec_2, tr_node * inserted_node){
	int hash_val, hash_val_mod;
	
	tr_node * ptr;
	//~ static int clashes = 0;
	//~ static int total = 0;
	
	/* where to insert */
	inserted_node->hash_value = hash_fci(type, loop_spec_1, loop_spec_2, inserted_node->seqlen);
	hash_val_mod = pow_of_2_modulus(inserted_node->hash_value, hash_t->hash_size);
	
	//~ total++;
	
	if(!hash_t->array[hash_val_mod]){ /* head case */
		hash_t->array[hash_val_mod] = inserted_node;
		return;
	}//else{
		//~ clashes++;	
	//~ }
	ptr = hash_t->array[hash_val_mod];
	while(ptr->next_in_hash){ /* inside list */
		ptr = ptr->next_in_hash;
	}
	//printf("Detected %d clashes out of %d cases for total %f ratio of cases.\n", clashes, total, clashes*1./total);
	ptr->next_in_hash = inserted_node;
}


/* creates and inserts a new tr_node (along with creating hash) into tree */
tr_node* insert_tr_node(int type, int loop_spec_1, int loop_spec_2, tr_node * last_node){  
  /* forge new tr_node */	  
  tr_node * new_tr_node = create_tr_node(type, loop_spec_1, loop_spec_2, last_node->seqlen, last_node);
  if (!last_node->child){ /* empty so no child yet - insert it directly */
	  last_node->child = new_tr_node;
  }else{
	  if(!last_node->hash_tab){ /* hash not defined yet (we have exactly two values) */
		  last_node->hash_tab = init_hash(last_node->seqlen);
		  insert_val(last_node->hash_tab, last_node->child->type, last_node->child->loop_spec_1, last_node->child->loop_spec_2, last_node->child);
		  last_node->hash_tab->hash_elemts++;
	  }
	  /* insert it into hashtable of previous node */
	  insert_val(last_node->hash_tab, type, loop_spec_1, loop_spec_2, new_tr_node);
	  last_node->hash_tab->hash_elemts++;
	  if(((double)last_node->hash_tab->hash_elemts/last_node->hash_tab->hash_size) > 0.5){
		  extend_hash(last_node->hash_tab);
	  }
  }
  return new_tr_node;
}


/* tracks whether there is node in hash with given values */
tr_node* access_val(tr_node *parent, int type, int loop_spec_1, int loop_spec_2){
	tr_node * ptr = NULL;
	if(!parent->hash_tab){ /* hash does not exist (at most one child) */
		if(parent->child){ /* child exists (exactly one child) */
			if((parent->child->type == type) && (parent->child->loop_spec_1 == loop_spec_1)
			 && (parent->child->loop_spec_2 == loop_spec_2)){
				return parent->child;
			}
		}
	}else{ /* hash exists -> 2+ children */
		int hash_f = hash_fci(type, loop_spec_1, loop_spec_2, parent->seqlen);
		int hash_f_mod = pow_of_2_modulus(hash_f, parent->hash_tab->hash_size);
		ptr = parent->hash_tab->array[hash_f_mod];
		while(ptr){
			if(ptr->hash_value == hash_f){	
				break;
			}
			ptr = ptr->next_in_hash;
		}
	}
	if(ptr){
		return ptr;	
	}
	return NULL;
}


/* creates node with properties that don't already exist for given node */
tr_node * add_if_nexists(int type, int loop_spec_1, int loop_spec_2, tr_node *last_node){
  tr_node * new_tr_node = access_val(last_node, type, loop_spec_1, loop_spec_2);
  if(!new_tr_node){
	new_tr_node = insert_tr_node(type, loop_spec_1, loop_spec_2, last_node);
  }
  return new_tr_node;  	
}


/* returns weight of node if it exists; otherwise returns 0 */
double tr_node_weight(tr_node * last_node, int type, int loop_spec_1, int loop_spec_2){
  tr_node * next_node = access_val(last_node, type, loop_spec_1, loop_spec_2);
  if(next_node)
	 return next_node->weight;  
 
  return 0.;
}


/* returns sum of all weights for given parent */ 
double total_weight_par(tr_node * last_node){
	if((!last_node->hash_tab)&&(!last_node->child)){
		return 0.;	
	}
	return last_node->weight;
}


/* returns sum of all weights for given parent and for given type*/
double total_weight_par_type(int type, tr_node * last_node){
  tr_node *ptr = NULL;
  double weight_sum = 0.;
  int i;
  if(!last_node->hash_tab){
	if(last_node->child){
		return last_node->child->weight;
	}
  }
  else{
	  for(i=0; i<last_node->hash_tab->hash_size; i++){
		 ptr = last_node->hash_tab->array[i];
		 while(ptr){
		   if(ptr->type == type){	   	 
			 weight_sum += ptr->weight;
		   } 
		   ptr = ptr->next_in_hash; 	
		 }
	  }
	  return weight_sum;
  }
}


/* updates weight of a given node */
void update_weight(tr_node * node, double weight){
  node->weight += weight;
}



/* tracebacks to root while updating values for each node passed through */
tr_node * traceback_to_root(tr_node * leaf, double struct_weight){
  update_weight(leaf, struct_weight);
  while(leaf->parent){  
	update_weight(leaf->parent, struct_weight);
    leaf = leaf->parent;
  }
  return leaf; 	
}


/* destructor */
void free_all_nr(tr_node * root){
	if(!root->hash_tab){
		if(root->child)
			free_all_nr(root->child);	
	}else{
		for(int i=0; i<root->hash_tab->hash_size; i++){
			tr_node * ptr =  root->hash_tab->array[i];
			tr_node * tmp; 
			while(ptr){
				tmp = ptr->next_in_hash;
				free_all_nr(ptr);
				ptr = tmp;
			}
		}
		free(root->hash_tab->array);
		free(root->hash_tab);
	}
	free(root);
}

/*********************************************************/
/*       data structure (linked_list + hash table)       */
/*********************************************************/

/* This creates structure that uses linked list instead of hash. The thought behind this is
 * the order of investigated nodes is always the same so we can add them to specific place.
 * It is thus a bit faster.
 */

tllr_node * create_tllr_node(int type, int loop_spec_1, int loop_spec_2, tllr_node * parent){
	static int memory_index = 0;	
	tllr_node * new_tllr_node = &nr_memory_allocated[memory_index];
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
	new_tllr_node->loop_spec_1 = loop_spec_1;
	new_tllr_node->loop_spec_2 = loop_spec_2;
	new_tllr_node->parent = parent;
	new_tllr_node->next_node = NULL;
	new_tllr_node->head = NULL;
	new_tllr_node->weight = 0;
	
	memory_index++;
	return new_tllr_node;   
}

/* creates root (start of a tree) */
tllr_node* create_ll_root(){
	tllr_node* root = create_tllr_node(NRT_NONE_TYPE, 0, 0, NULL);
	return root;  	
}

/* inserts a tllr_node before 'next_node' and after previous node
 * Node cursor hols previous and current ll_node */
tllr_node* insert_tllr_node(tllr_node* memorized_node_prev, tllr_node* memorized_node_cur, 
							int type, int loop_spec_1, int loop_spec_2, tllr_node* parent_node){
	tllr_node * new_node = create_tllr_node(type, loop_spec_1, loop_spec_2, parent_node);
	
	if(!memorized_node_prev){ /* first node to be inserted */
		parent_node->head = new_node;	
	}else{
		memorized_node_prev->next_node = new_node;
	}
	new_node->next_node = memorized_node_cur;
	return new_node;
}

/* resets cursor to beginning of loop*/
void reset_cursor(tllr_node** memorized_node_prev, tllr_node** memorized_node_cur, tllr_node* current_node){
	(*memorized_node_prev) = NULL;
	(*memorized_node_cur) = current_node->head;
}

/* advances pointer in loop if the identifier coincide with current pointer and returns weight */
inline void advance_cursor(tllr_node** memorized_node_prev, tllr_node** memorized_node_cur, int type, int loop_spec_1, int loop_spec_2){
	if(*memorized_node_cur){
		if((*memorized_node_cur)->type == type 
		 && (*memorized_node_cur)->loop_spec_1 == loop_spec_1 
		 && (*memorized_node_cur)->loop_spec_2 == loop_spec_2){

			(*memorized_node_prev) = (*memorized_node_cur);
			(*memorized_node_cur) = (*memorized_node_cur)->next_node;
		}
	}
}

/* gets weight of actual node */
inline double get_weight(tllr_node* memorized_node_cur, int type, int loop_spec_1, int loop_spec_2){
	double weight = 0;
	if(memorized_node_cur){
		if(memorized_node_cur->type == type 
			 && memorized_node_cur->loop_spec_1 == loop_spec_1 
			 && memorized_node_cur->loop_spec_2 == loop_spec_2){
				 weight = memorized_node_cur->weight;		
		}
	}
	return weight;
}

/* get weight of all child nodes */
double get_weight_all(tllr_node* last_node){
	if(!last_node->head){
		return 0;	
	}	
	return last_node->weight;
}

/* get weight of all child nodes of certain type */
double get_weight_type_spec(int type, tllr_node* last_node){
	double weight_total = 0;
	tllr_node* ptr = last_node->head;
	while(ptr){
		if(ptr->type == type){
			weight_total += ptr->weight;
		}
		ptr = ptr->next_node;
	}
	return weight_total;
}

/* adds node if the current one isn't the one we want */
inline tllr_node* add_if_nexists_ll(int type, int loop_spec_1, int loop_spec_2, tllr_node* memorized_node_prev, tllr_node* memorized_node_cur, tllr_node* parent_node){
	tllr_node* returned_node;
	if(memorized_node_cur){
		if(memorized_node_cur->type == type 
		 && memorized_node_cur->loop_spec_1 == loop_spec_1 
		 && memorized_node_cur->loop_spec_2 == loop_spec_2){
			 returned_node = memorized_node_cur;
		}else{
			 returned_node = insert_tllr_node(memorized_node_prev, memorized_node_cur, type, loop_spec_1, loop_spec_2, parent_node);
		}
	}else{
		 returned_node = insert_tllr_node(memorized_node_prev, memorized_node_cur, type, loop_spec_1, loop_spec_2, parent_node);
	}
	return returned_node;
}

/* updates weight of a given node */
void update_weight_ll(tllr_node * node, double weight){
  node->weight += weight;
}

/* tracebacks to root while updating values for each node passed through */
tllr_node * traceback_to_ll_root(tllr_node * leaf, double weight){
  update_weight_ll(leaf, weight);
  while(leaf->parent){  
	update_weight_ll(leaf->parent, weight);
    leaf = leaf->parent;
  }
  return leaf; 	
}

/* destructor */
void free_all_ll(tllr_node * root){
	tllr_node* ptr;
	if(root->head){
		free_all_ll(root->head);
		if(root->next_node){
			 free_all_ll(root->next_node);	
		}
		free(root);	
	} 
}