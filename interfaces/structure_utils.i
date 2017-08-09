/*#######################################*/
/* Structure Utitlities section          */
/*#######################################*/

/* compressing / decompressing dot-bracket strings */
%rename (db_pack) vrna_db_pack;
%newobject vrna_db_pack;

%ignore pack_structure;
%rename (pack_structure) my_pack_structure;
%newobject my_pack_structure;
%{
  char *my_pack_structure(const char *s){
    return vrna_db_pack(s);
  }
%}
char *my_pack_structure(const char *s);

%rename (db_unpack) vrna_db_unpack;
%newobject vrna_db_unpack;

%ignore unpack_structure;
%rename (unpack_structure) my_unpack_structure;
%newobject my_unpack_structure;
%{
  char *my_unpack_structure(const char *packed){
    return vrna_db_unpack(packed);
  }
%}
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

%{
#include <vector>

  std::vector<int> my_ptable(std::string str){
    short int* pt = vrna_ptable(str.c_str());
    std::vector<int> v_pt;
    int i;

    for(i=0; i <= pt[0]; i++){
      v_pt.push_back(pt[i]);
    }
    free(pt);
    return v_pt;
  }
%}

std::vector<int> my_ptable(std::string str);


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

char *my_db_from_ptable(std::vector<int> pt);


/* pair table related functions */
%ignore make_pair_table;
%ignore make_pair_table_pk;
%ignore copy_pair_table;
%ignore alimake_pair_table;
%ignore make_pair_table_snoop;
%ignore make_loop_index_pt;

%ignore vrna_hx_t;
%ignore vrna_hx_s;

%ignore bp_distance;
%rename (bp_distance) my_bp_distance;
%{
  int my_bp_distance(const char *str1, const char *str2){
    return vrna_bp_distance(str1,str2);
  }
%}
int my_bp_distance(const char *str1, const char *str2);

%ignore make_referenceBP_array;
%ignore compute_BPdifferences;
%ignore parenthesis_structure;
%ignore parenthesis_zuker;
%ignore letter_structure;
%ignore bppm_to_structure;
%ignore bppm_symbol;


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
  %extend {
    char *__str__() {
      char *tmp = vrna_strdup_printf("[ i: %d, j: %d, p: %.10g, t: %2d ]", $self->i, $self->j, $self->p, $self->type);
      return tmp;
    }
  }
} vrna_ep_t;

char *vrna_ep_t::__str__();


/*
 * wrap functions that deal with 'vrna_ep_t' I/O
 */
%ignore plist;
%ignore assign_plist_from_db;
%ignore assign_plist_from_pr;

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

std::vector<vrna_ep_t> my_plist(std::string structure, float pr);
std::string db_from_plist(std::vector<vrna_ep_t> elem_probs, unsigned int length);


%extend vrna_fold_compound_t {
#include <vector>
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

%include  <ViennaRNA/structure_utils.h>
