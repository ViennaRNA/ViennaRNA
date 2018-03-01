#ifndef DATA_STRUCTURES_NONRED_H_
#define DATA_STRUCTURES_NONRED_H_

/* explanation: each node (tr_node) holds a pointer to a start of a linked
 * list to the next node, link to previous node, Boltzmann's factor
 * for given flat structure and itentifier, which corresponds to an id of
 * chosen flat structure
 */

/* define node of a tree */
typedef struct tr_node tr_node; 

/* defines types for different sequence stetches/loops */
typedef enum _type Type;
enum stetch_type{
	    NRT_NONE_TYPE = 0, 			/* nonetype: both 0 (root) */
	    NRT_HAIRPIN = 1,			/* hairpin: both 0 (unused) */
	    NRT_IT_LOOP = 2,			/* internal loop: (i', j') - closed pair */
	    NRT_MT_LOOP = 3,			/* multiloop: K (second is 0 - unused) */
	    NRT_EXT_LOOP = 4,			/* external loop: (i, j) */
	    NRT_UNPAIRED_SG = 5,		/* unpaired positions: (first unpaired position - end known of stretch) */
	    NRT_QM1_BRANCH = 6,			/* branch in multiloop in qm1: (i,j) */
	    NRT_QM_PAIR = 7,			/* QM split on K with pairing within [i..k-1] : (k, 0)  */
	    NRT_QM_UNPAIR = 8,			/* QM split on K with no pairing within [i..k-1] : (k, 0) */
	};


/* define elements of hashtable */
typedef struct hash_node hash_node; /* list of arrays sorting triplets for every hashvalue*/
		
struct tr_node{
			int type; 
			int loop_spec_1; 
			int loop_spec_2; 
			int seqlen;
			tr_node * parent; 
			tr_node * child;  /* used when this node has at most one child */
			hash_node * hash_tab; /* used when there is at least two children */
			tr_node * next_in_hash; /* used for node that is next in hash table */
			int hash_value; /* hash value computed for this node */
			double weight;
		};

struct hash_node{
			tr_node **array;
			int hash_elemts;
			int hash_size;
		};
		
/* tree + hash functions */
/** @brief creates a root of datastructure tree (hash table version) **/
tr_node * create_root(int seqlen);

/** @brief returns a weigh of node (type, loop_spec_1, loop_spec_2) if child of last_node, otherwise returns 0.0 **/
double tr_node_weight(tr_node * last_node, int type, int loop_spec_1, int loop_spec_2);

/** @brief sums weight of all children of par_node and returns it **/
double total_weight_par(tr_node * par_node);

/** @brief sums weight of all children of par_node with certain type and returns it **/
double total_weight_par_type(int type, tr_node * par_node);

/** @brief creates node (type, loop_spec_1, loop_spec_2) if not existing and returns pointer to it, 
 * or returns pointer to exisiting case **/
tr_node * add_if_nexists(int type, int loop_spec_1, int loop_spec_2, tr_node*par_node);

/** @brief traces back from leaf to root while updating weights of leaf to all nodes in path, 
 *  returns pointer to root **/
tr_node * traceback_to_root(tr_node * leaf, double struct_weight);

/** @brief destructor **/
void free_all_nr(tr_node * root);

/******************************************/
/*       version with linked lists        */
/******************************************/   

typedef struct tllr_node tllr_node;

struct tllr_node{
			int type;
			int loop_spec_1;
			int loop_spec_2;
			tllr_node* parent;		/* vertical chaining - ancestor */
			tllr_node* head;		/* vertical chaining - successor */
			tllr_node* next_node;	/*horizontal chaining - linked list */
			double weight;
};

/* global variable for storing global nodes */
tllr_node *nr_memory_allocated; 

/* tree + linked list functions */
/** @brief creates a root of datastructure tree (linked list version) **/
tllr_node * create_ll_root();

/** resets cursor to current_node and start of linked list **/
void reset_cursor(tllr_node** memorized_node_prev, tllr_node** memorized_node_cur, tllr_node* current_node);

/** @brief moves cursor to next node if current_node is identical to one in loop, otherwise does nothing **/
void advance_cursor(tllr_node** memorized_node_prev, tllr_node** memorized_node_cur, int type, int loop_spec_1, int loop_spec_2);

/** @brief returns a weigh of node (type, loop_spec_1, loop_spec_2) if child of last_node, otherwise returns 0.0 **/
double get_weight(tllr_node* memorized_node_cur, int type, int loop_spec_1, int loop_spec_2);

/** @brief sums weight of all children of par_node and returns it **/
double get_weight_all(tllr_node * par_node);

/** @brief sums weight of all children of par_node with certain type and returns it **/
double get_weight_type_spec(int type, tllr_node * par_node);

/** @brief creates node (type, loop_spec_1, loop_spec_2) if not existing and returns pointer to it, 
 * or returns pointer to exisiting case **/
tllr_node* add_if_nexists_ll(int type, int loop_spec_1, int loop_spec_2, tllr_node* memorized_node_prev, tllr_node* memorized_node_cur, tllr_node* parent_node);

/** @brief traces back from leaf to root while updating weights of leaf to all nodes in path, 
 *  returns pointer to root **/
tllr_node * traceback_to_ll_root(tllr_node * leaf, double weight);

/** @brief destructor **/
void free_all_ll(tllr_node * root);


#endif