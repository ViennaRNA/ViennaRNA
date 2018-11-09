/*#######################################*/
/* Structure Utitlities section          */
/*#######################################*/

/* first, the ignore list */
%ignore pack_structure;
%ignore unpack_structure;
%ignore make_pair_table;
%ignore make_pair_table_pk;
%ignore copy_pair_table;
%ignore alimake_pair_table;
%ignore make_pair_table_snoop;
%ignore make_loop_index_pt;
%ignore bp_distance;
%ignore make_referenceBP_array;
%ignore compute_BPdifferences;
%ignore parenthesis_structure;
%ignore parenthesis_zuker;
%ignore letter_structure;
%ignore bppm_to_structure;
%ignore bppm_symbol;
%ignore plist;
%ignore assign_plist_from_db;
%ignore assign_plist_from_pr;

%ignore vrna_hx_t;
%ignore vrna_hx_s;
%ignore vrna_ptable_from_string;
%ignore vrna_db_flatten;
%ignore vrna_db_flatten_to;
%ignore vrna_db_from_WUSS;


/* next, all the wrappers */


/* compressing / decompressing dot-bracket strings */
%rename (db_pack) vrna_db_pack;
%newobject vrna_db_pack;

%rename (pack_structure) my_pack_structure;
%newobject my_pack_structure;

%{
  char *my_pack_structure(const char *s){
    return vrna_db_pack(s);
  }
%}
#ifdef SWIGPYTHON
%feature("autodoc") my_pack_structure;
%feature("kwargs") my_pack_structure;
#endif

char *my_pack_structure(const char *s);

%rename (db_unpack) vrna_db_unpack;
%newobject vrna_db_unpack;

%rename (unpack_structure) my_unpack_structure;
%newobject my_unpack_structure;
%{
  char *my_unpack_structure(const char *packed){
    return vrna_db_unpack(packed);
  }
%}
#ifdef SWIGPYTHON
%feature("autodoc") my_unpack_structure;
%feature("kwargs") my_unpack_structure;
#endif

char *my_unpack_structure(const char *packed);

%inline %{
  short *make_loop_index(const char *structure) {
  /* number each position by which loop it belongs to (positions start at 0) */
    int i,hx,l,nl;
    int length;
    short *stack;
    short *loop;
    length = strlen(structure);
    stack = (short *) vrna_alloc(sizeof(short)*(length+1));
    loop = (short *) vrna_alloc(sizeof(short)*(length+2));
    hx=l=nl=0;
    for (i=0; i<length; i++) {
      if (structure[i] == '(') {
        nl++; l=nl;
        stack[hx++]=i;
      }
      loop[i]=l;
      if (structure[i] ==')') {
        --hx;
        if (hx>0)
          l = loop[stack[hx-1]];  /* index of enclosing loop   */
        else l=0;                 /* external loop has index 0 */
        if (hx<0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in make_loop_index");
        }
      }
    }
    free(stack);
    return loop;
  }
%}

%rename (ptable) my_ptable;
%rename (ptable_from_string) my_ptable_from_string;

%{
#include <vector>

  std::vector<int>
  my_ptable(std::string str)
  {
    short int* pt = vrna_ptable(str.c_str());
    std::vector<int> v_pt;
    int i;

    for(i=0; i <= pt[0]; i++){
      v_pt.push_back(pt[i]);
    }
    free(pt);
    return v_pt;
  }

  std::vector<int>
  my_ptable_from_string(std::string str,
                        unsigned int options = VRNA_BRACKETS_DEFAULT)
  {
    short int         *pt;
    int               i;
    std::vector<int>  v_pt;

    pt = vrna_ptable_from_string(str.c_str(), options);

    for(i = 0; i <= pt[0]; i++)
      v_pt.push_back(pt[i]);

    free(pt);
    return v_pt;
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_ptable;
%feature("kwargs") my_ptable;
%feature("autodoc") my_ptable_from_string;
%feature("kwargs") my_ptable_from_string;
#endif

std::vector<int> my_ptable(std::string str);
std::vector<int> my_ptable_from_string(std::string str, unsigned int options = VRNA_BRACKETS_DEFAULT);

/* pairtable with pseudoknots */
%rename (ptable_pk) my_ptable_pk;

%{
#include <vector>

  std::vector<int> my_ptable_pk(std::string str){
   short int* pt_pk = vrna_pt_pk_get(str.c_str());
    std::vector<int> v_pt;
    int i;

    for(i=0; i <= pt_pk[0]; i++){
      v_pt.push_back(pt_pk[i]);
    }
    free(pt_pk);
    return v_pt; 
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_ptable_pk;
%feature("kwargs") my_ptable_pk;
#endif

std::vector<int> my_ptable_pk(std::string str);


%rename (db_from_ptable) my_db_from_ptable;

%{
#include <vector>
  char *my_db_from_ptable(std::vector<int> pt)
  {
    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
    return vrna_db_from_ptable((short*)&vc[0]);
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_db_from_ptable;
%feature("kwargs") my_db_from_ptable;
#endif

char *my_db_from_ptable(std::vector<int> pt);


%rename (loopidx_from_ptable) my_loopidx_from_ptable;

%{

  std::vector<int>
  my_loopidx_from_ptable(std::vector<int> pt)
  {
    int                 i, *idx;
    std::vector<short>  vc;
    std::vector<int>    v_idx;

    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);

    idx = vrna_loopidx_from_ptable((short *)&vc[0]);

    v_idx.assign(idx, idx + pt.size());

    free(idx);

    return v_idx;
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_loopidx_from_ptable;
%feature("kwargs") my_loopidx_from_ptable;
#endif

std::vector<int> my_loopidx_from_ptable(std::vector<int> pt);


%rename (bp_distance) my_bp_distance;
%{
  int my_bp_distance(const char *str1, const char *str2){
    return vrna_bp_distance(str1,str2);
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_bp_distance;
%feature("kwargs") my_bp_distance;
#endif

int my_bp_distance(const char *str1, const char *str2);

%rename (dist_mountain) my_dist_mountain;

%{
  double
  my_dist_mountain( std::string str1,
                    std::string str2,
                    unsigned int p = 1)
  {
    return vrna_dist_mountain(str1.c_str(), str2.c_str(), p);
  }

%}

double
my_dist_mountain(std::string str1,
                 std::string str2,
                 unsigned int p = 1);


/*
 *  Rename 'struct vrna_ep_t' to target language class 'ep'
 *  and create necessary wrapper extensions
 */
%rename (ep) vrna_ep_t;

%newobject vrna_ep_t::__str__;

typedef struct {
  int i;
  int j;
  float p;
  int type;
} vrna_ep_t;

%extend vrna_ep_t {
    char *__str__() {
      char *tmp = vrna_strdup_printf("[ i: %d, j: %d, p: %.10g, t: %2d ]", $self->i, $self->j, $self->p, $self->type);
      return tmp;
    }
}


/*
 * wrap functions that deal with 'vrna_ep_t' I/O
 */

%rename (plist) my_plist;

%{
#include <vector>
  std::vector<vrna_ep_t> my_plist(std::string structure, float pr = 0.95*0.95)
  {
    std::vector<vrna_ep_t > ep_v;
    vrna_ep_t               *ptr, *plist;

    plist = vrna_plist(structure.c_str(), pr);

    for (ptr = plist; ptr->i && ptr->j; ptr++) {
      vrna_ep_t pl;
      pl.i = ptr->i;
      pl.j = ptr->j;
      pl.p = ptr->p;
      pl.type = ptr->type;
      ep_v.push_back(pl);
    }

    free(plist);

    return ep_v;
  }

  std::string db_from_plist(std::vector<vrna_ep_t> pairs, unsigned int length) {
    char *str = vrna_db_from_plist(&pairs[0], length);
    std::string ret(str);
    free(str);
    return ret;
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_plist;
%feature("kwargs") my_plist;
%feature("autodoc") db_from_plist;
%feature("kwargs") db_from_plist;
#endif

std::vector<vrna_ep_t> my_plist(std::string structure, float pr);
std::string db_from_plist(std::vector<vrna_ep_t> elem_probs, unsigned int length);


%extend vrna_fold_compound_t {
#include <vector>

#ifdef SWIGPYTHON
%feature("autodoc") plist_from_probs;
%feature("kwargs") plist_from_probs;
#endif

  std::vector<vrna_ep_t> plist_from_probs(double cutoff) {
    std::vector<vrna_ep_t > ep_v;
    vrna_ep_t               *ptr, *plist;

    plist = vrna_plist_from_probs($self, cutoff);

    for (ptr = plist; ptr->i && ptr->j; ptr++) {
      vrna_ep_t pl;
      pl.i = ptr->i;
      pl.j = ptr->j;
      pl.p = ptr->p;
      pl.type = ptr->type;
      ep_v.push_back(pl);
    }

    free(plist);

    return ep_v;
  }
}


/*
 *  Wrap constants required for 'vrna_ep_t' interaction
 */
%constant int PLIST_TYPE_BASEPAIR = VRNA_PLIST_TYPE_BASEPAIR;
%constant int PLIST_TYPE_GQUAD    = VRNA_PLIST_TYPE_GQUAD;
%constant int PLIST_TYPE_H_MOTIF  = VRNA_PLIST_TYPE_H_MOTIF;
%constant int PLIST_TYPE_I_MOTIF  = VRNA_PLIST_TYPE_I_MOTIF;
%constant int PLIST_TYPE_UD_MOTIF = VRNA_PLIST_TYPE_UD_MOTIF;
%constant int PLIST_TYPE_STACK    = VRNA_PLIST_TYPE_STACK;


%constant unsigned int STRUCTURE_TREE_HIT            = VRNA_STRUCTURE_TREE_HIT;
%constant unsigned int STRUCTURE_TREE_SHAPIRO_SHORT  = VRNA_STRUCTURE_TREE_SHAPIRO_SHORT;
%constant unsigned int STRUCTURE_TREE_SHAPIRO        = VRNA_STRUCTURE_TREE_SHAPIRO;
%constant unsigned int STRUCTURE_TREE_SHAPIRO_EXT    = VRNA_STRUCTURE_TREE_SHAPIRO_EXT;
%constant unsigned int STRUCTURE_TREE_SHAPIRO_WEIGHT = VRNA_STRUCTURE_TREE_SHAPIRO_WEIGHT;
%constant unsigned int STRUCTURE_TREE_EXPANDED       = VRNA_STRUCTURE_TREE_EXPANDED;


/*
 *  Wrap more library routines for secondary structure representation formats
 *  dot-bracket, WUSS, etc.
 */
%{

  void db_flatten(char *structure, unsigned int options = VRNA_BRACKETS_DEFAULT)
  {
    vrna_db_flatten(structure, options);
  }

  void db_flatten(char *structure, std::string target, unsigned int options = VRNA_BRACKETS_DEFAULT)
  {
    if (target.size() == 2)
      vrna_db_flatten_to(structure, target.c_str(), options);
    else
      vrna_message_warning("db_flatten(): target pair must be string of exactly 2 characters!");
  }

  std::string
  db_from_WUSS(std::string wuss)
  {
    char *c_str = vrna_db_from_WUSS(wuss.c_str());
    std::string db = c_str;
    free(c_str);
    return db;
  }

  std::string
  db_to_tree_string(std::string   structure,
                    unsigned int  type)
  {
    char *c_str = vrna_db_to_tree_string(structure.c_str(), type);
    std::string tree = c_str;
    free(c_str);
    return tree;
  }

  std::string
  tree_string_unweight(std::string structure)
  {
    char *c_str = vrna_tree_string_unweight(structure.c_str());
    std::string tree = c_str;
    free(c_str);
    return tree;
  }

  std::string
  tree_string_to_db(std::string structure)
  {
    char *c_str = vrna_tree_string_to_db(structure.c_str());
    std::string db = c_str;
    free(c_str);
    return db;
  }

%}

#ifdef SWIGPYTHON
%feature("autodoc") db_to_tree_string;
%feature("kwargs") db_to_tree_string;
#endif

void        db_flatten(char *structure, unsigned int options = VRNA_BRACKETS_DEFAULT);
void        db_flatten(char *structure, std::string target, unsigned int options = VRNA_BRACKETS_DEFAULT);
std::string db_from_WUSS(std::string wuss);
std::string db_to_tree_string(std::string structure, unsigned int type);
std::string tree_string_unweight(std::string structure);
std::string tree_string_to_db(std::string structure);


%constant unsigned int BRACKETS_RND   = VRNA_BRACKETS_RND;
%constant unsigned int BRACKETS_ANG   = VRNA_BRACKETS_ANG;
%constant unsigned int BRACKETS_SQR   = VRNA_BRACKETS_SQR;
%constant unsigned int BRACKETS_CLY   = VRNA_BRACKETS_CLY;
%constant unsigned int BRACKETS_ALPHA = VRNA_BRACKETS_ALPHA;
%constant unsigned int BRACKETS_DEFAULT = VRNA_BRACKETS_DEFAULT;

%newobject vrna_db_to_element_string;
%rename (db_to_element_string) vrna_db_to_element_string;


%include  <ViennaRNA/utils/structures.h>
