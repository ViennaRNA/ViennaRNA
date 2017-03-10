/*#######################################*/
/* Utitlities section                    */
/*#######################################*/

/**********************************************/
/* BEGIN interface for generic utilities      */
/**********************************************/

%ignore get_line;
%ignore get_input_line;
%ignore get_ptypes;
%ignore get_indx;
%ignore get_iindx;
%ignore print_tty_input_seq;
%ignore print_tty_input_seq_str;
%ignore warn_user;
%ignore nrerror;
%ignore space;
%ignore xrealloc;
%ignore init_rand;
%ignore urn;
%ignore int_urn;
%ignore filecopy;
%ignore time_stamp;


%include  <ViennaRNA/utils.h>

/**********************************************/
/* BEGIN interface for string utilities       */
/**********************************************/

/* random string */
%ignore random_string;
%rename (random_string) vrna_random_string;
%newobject vrna_random_string;

/* hamming distance */
%rename (hamming_distance) vrna_hamming_distance;
%rename (hamming_distance_bound) vrna_hamming_distance_bound;

%ignore hamming;
%ignore hamming_bound;

%rename (hamming) my_hamming;
%{
  int my_hamming(const char *s1, const char *s2){
    return vrna_hamming_distance(s1, s2);
  }
%}
int my_hamming(const char *s1, const char *s2);

%rename (hamming_bound) my_hamming_bound;
%{
  int my_hamming_bound(const char *s1, const char *s2, int n){
    return vrna_hamming_distance_bound(s1, s2, n);
  }
%}
int my_hamming_bound(const char *s1, const char *s2, int n);

/* RNA -> DNA conversion */
%ignore str_DNA2RNA;

/* string uppercase
 * (there is surely a more efficient version in the target language,
 * so we do not wrap them)
 */
%ignore str_uppercase;

/* encoding / decoding of nucleotide sequences */

%{

#include <cstring>

short *encode_seq(char *sequence) {
  unsigned int i,l;
  short *S;
  l = strlen(sequence);
  S = (short *) vrna_alloc(sizeof(short)*(l+2));
  S[0] = (short) l;

  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++)
    S[i]= (short) encode_char(toupper(sequence[i-1]));

  /* for circular folding add first base at position n+1 */
  S[l+1] = S[1];

  return S;
}
%}
short *encode_seq(char *sequence);

%include  <ViennaRNA/string_utils.h>

/**********************************************/
/* BEGIN interface for structure utilities    */
/**********************************************/

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

%ignore assign_plist_from_db;
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
%ignore assign_plist_from_pr;
%ignore parenthesis_structure;
%ignore parenthesis_zuker;
%ignore letter_structure;
%ignore bppm_to_structure;
%ignore bppm_symbol;


typedef struct {
  int i;
  int j;
  float p;
  int type;
} vrna_plist_t;

%nodefaultdtor vrna_plist_t;

%extend vrna_plist_t {

  ~vrna_plist_t() {
    free($self);
  }
}

namespace std {
  %template(PlistVector) std::vector<vrna_plist_t>;
};

%rename (plist) my_plist;

%{
#include <vector>
  std::vector<vrna_plist_t> my_plist(std::string structure, float pr = 0.95*0.95){
    
    std::vector<vrna_plist_t > vc;
    vrna_plist_t *ptr, *plist;
    
    plist = vrna_plist(structure.c_str(),pr);

    for(ptr=plist; ptr->i && ptr->j; ptr++){
      vrna_plist_t pl;
      pl.i = ptr->i;
      pl.j = ptr->j;
      pl.p = ptr->p;
      pl.type = ptr->type;
      vc.push_back(pl);
    }

    free(plist);

    return vc;
  }
%}

%newobject my_plist;
std::vector<vrna_plist_t> my_plist(std::string structure, float pr);

%constant int PLIST_TYPE_BASEPAIR = VRNA_PLIST_TYPE_BASEPAIR;
%constant int PLIST_TYPE_GQUAD    = VRNA_PLIST_TYPE_GQUAD;
%constant int PLIST_TYPE_H_MOTIF  = VRNA_PLIST_TYPE_H_MOTIF;
%constant int PLIST_TYPE_I_MOTIF  = VRNA_PLIST_TYPE_I_MOTIF;
%constant int PLIST_TYPE_UD_MOTIF = VRNA_PLIST_TYPE_UD_MOTIF;

%include  <ViennaRNA/structure_utils.h>

/**********************************************/
/* BEGIN interface for alignment utilities    */
/**********************************************/


%rename (aln_consensus) my_consensus;
%rename (consensus) my_consensus;
%{
#include <vector>

  std::string my_consensus(std::vector<std::string> alignment){

    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    char *c = consensus((const char **)&v[0]);
    std::string cons(c);
    free(c);
    return cons;
  }

%}

std::string my_consensus(std::vector<std::string> alignment);

%ignore consensus;


%rename (aln_consensus_mis) my_consensus_mis;
%rename (consens_mis) my_consensus_mis;
%{
#include <vector>

  std::string my_consensus_mis(std::vector<std::string> alignment){

    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    char *c = consens_mis((const char **)&v[0]);
    std::string cons(c);
    free(c);
    return cons;
  }

%}

std::string my_consensus_mis(std::vector<std::string> alignment);

%ignore consens_mis;


%rename (aln_mpi) my_aln_mpi;
%{
#include <vector>

  int my_aln_mpi(std::vector<std::string> alignment){

    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    int mpi = vrna_aln_mpi((const char **)&v[0]);

    return mpi;
  }

%}

int my_aln_mpi(std::vector<std::string> alignment);


%rename (aln_pscore) my_aln_pscore;

%{
#include <vector>

  std::vector<std::vector<int> > my_aln_pscore(std::vector<std::string> alignment, vrna_md_t *md = NULL){

    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    std::vector<std::vector<int> > pscore;
    int *ps = vrna_aln_pscore((const char **)&v[0], md);

    int n     = alignment[0].length();
    int *idx  = vrna_idx_col_wise(n);

    std::vector<int> z_row(n+1, 0);
    pscore.push_back(z_row);

    for(int i = 1; i < n; i++){
      std::vector<int> score_i;
      score_i.push_back(0);
      for(int j = 1; j <= i; j++)
        score_i.push_back(ps[idx[i] + j]);
      for(int j = i + 1; j <= n; j++)
        score_i.push_back(ps[idx[j] + i]);
      pscore.push_back(score_i);
    }

    free(ps);
    free(idx);

    return pscore;
  }

%}

std::vector<std::vector<int> > my_aln_pscore(std::vector<std::string> alignment, vrna_md_t *md = NULL);


%ignore read_clustal;
%ignore get_ungapped_sequence;
%ignore get_mpi;
%ignore encode_ali_sequence;
%ignore alloc_sequence_arrays;
%ignore free_sequence_arrays;

%include  <ViennaRNA/aln_util.h>

/**********************************************/
/* BEGIN interface for Move_Set utilities    */
/**********************************************/

%rename (move_standard) my_move_standard;

%{
  char *my_move_standard(int *OUTPUT, char *seq, char *struc, enum MOVE_TYPE type,int verbosity_level, int shifts, int noLP){
    char *structure =  (char *)calloc(strlen(struc)+1,sizeof(char));
    strcpy(structure,struc);
    *OUTPUT = move_standard(seq,structure,type,verbosity_level,shifts,noLP);
    return structure;   
  }
%}
%newobject my_move_standard;
char *my_move_standard(int *OUTPUT, char *seq, char *struc, enum MOVE_TYPE type,int verbosity_level, int shifts, int noLP);
%ignore move_standard;

%include  <ViennaRNA/move_set.h>


/**********************************************/
/* BEGIN interface for File utilities         */
/**********************************************/

%rename (filename_sanitize) my_filename_sanitize;

%{
  std::string my_filename_sanitize(std::string name) {
    std::string s;
    char *name_sanitized = vrna_filename_sanitize(name.c_str(), NULL);
    if (name_sanitized)
      s = (const char *)name_sanitized;
    free(name_sanitized);
    return s;
  }

  std::string my_filename_sanitize(std::string name, char c) {
    std::string s;
    char *name_sanitized = vrna_filename_sanitize(name.c_str(), &c);
    if (name_sanitized)
      s = (const char *)name_sanitized;
    free(name_sanitized);
    return s;
  }
%}

std::string my_filename_sanitize(std::string name);
std::string my_filename_sanitize(std::string name, char c);

%include  <ViennaRNA/file_utils.h>



