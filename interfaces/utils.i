/*#######################################*/
/* Utitlities section                    */
/*#######################################*/

/**********************************************/
/* BEGIN interface for generic utilities      */
/**********************************************/

%newobject space;
%newobject time_stamp;
%newobject get_line;

%include  "../src/ViennaRNA/utils.h"

/**********************************************/
/* BEGIN interface for string utilities       */
/**********************************************/

/* random string */
%ignore random_string;
%ignore vrna_random_string;
%rename (random_string) vrna_random_string;
%newobject random_string;

/* hamming distance */
%ignore vrna_hamming_distance;
%ignore vrna_hamming_distance_bound;
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
%ignore vrna_seq_toRNA;
%ignore str_DNA2RNA;

/* string uppercase
 * (there is surely a more efficient version in the target language,
 * so we do not wrap them)
 */
%ignore vrna_seq_toupper;
%ignore str_uppercase;

/* encoding / decoding of nucleotide sequences */

%{
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

%ignore vrna_seq_encode;
%ignore vrna_seq_encode_simple;
%ignore vrna_nucleotide_encode;
%ignore vrna_nucleotide_decode;
%ignore vrna_aln_encode;

/* insertion/removal of cut point '&' character in string */
%ignore vrna_cut_point_insert;
%ignore vrna_cut_point_remove;

%include  "../src/ViennaRNA/string_utils.h"

/**********************************************/
/* BEGIN interface for structure utilities    */
/**********************************************/

/* compressing / decompressing dot-bracket strings */
%ignore vrna_db_pack;
%rename (db_pack) vrna_db_pack;
%newobject db_pack;

%ignore pack_structure;
%rename (pack_structure) my_pack_structure;
%newobject pack_structure;
%{
  char *my_pack_structure(const char *s){
    return vrna_db_pack(s);
  }
%}
char *my_pack_structure(const char *s);

%ignore vrn_db_unpack;
%rename (db_unpack) vrna_db_unpack;
%newobject db_unpack;

%ignore unpack_structure;
%rename (unpack_structure) my_unpack_structure;
%newobject unpack_structure;
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

/* pair table related functions */
%ignore vrna_ptable;
%ignore vrna_pt_pk_get;
%ignore vrna_ptable_copy;
%ignore vrna_pt_ali_get;
%ignore vrna_pt_snoop_get;
%ignore vrna_db_from_ptable;
%ignore vrna_loopidx_from_ptable;
%ignore make_pair_table;
%ignore make_pair_table_pk;
%ignore copy_pair_table;
%ignore alimake_pair_table;
%ignore make_pair_table_snoop;
%ignore make_loop_index_pt;



%include  "../src/ViennaRNA/structure_utils.h"

/**********************************************/
/* BEGIN interface for alignment utilities    */
/**********************************************/


%newobject consensus;
%newobject consensus_mis;
char *consensus(const char **AS);
char *consens_mis(const char **AS);

