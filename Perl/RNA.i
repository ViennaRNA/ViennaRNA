%module RNA
//%pragma(perl5)  modulecode="@EXPORT=qw(fold);"
%pragma(perl5)  include="RNA.pod"
%{
#include  "../H/utils.h"
#include  "../H/fold_vars.h"
#undef fold
#include  "../H/fold.h"
#include  "../H/cofold.h"
#include  "../H/part_func.h"
#include  "../H/part_func_co.h"
#include  "../H/PS_dot.h"
#include  "../H/inverse.h"
#include  "../H/RNAstruct.h"
#include  "../H/treedist.h"
#include  "../H/stringdist.h"
#include  "../H/profiledist.h"
#include  "../H/dist_vars.h"
#include  "../H/pair_mat.h"
#include  "../H/subopt.h"
#include  "../H/energy_const.h"
#include  "../H/params.h"
#include  "../H/duplex.h"
#include  "../H/alifold.h"
#include  "../H/aln_util.h"
extern char *pbacktrack(char *seq);
extern void  read_parameter_file(const char fname[]);
extern void  write_parameter_file(const char fname[]);
extern int simple_xy_coordinates(short *pair_table, float *X, float *Y);
extern int naview_xy_coordinates(short *pair_table, float *X, float *Y);
%}
//
%include carrays.i
%array_functions(int, intP);
%array_class(int, intArray);
%array_functions(float, floatP);
%array_class(float, floatArray);
%array_functions(double, doubleP);
%array_class(double, doubleArray);
%array_functions(unsigned short, shortP);
%include cdata.i

%constant double VERSION = 0.3;
%include typemaps.i
%include tmaps.i  // additional typemaps

//%title "Interface to the Vienna RNA library"
//%section "Folding Routines"
//%subsection "Minimum free Energy Folding"

%rename (fold) my_fold;

%{
  char *my_fold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = fold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_fold;
char *my_fold(char *string, char *constraints = NULL, float *OUTPUT);
%ignore fold;
%include  "../H/fold.h"

%rename (cofold) my_cofold;

%{
  char *my_cofold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = cofold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_cofold;
char *my_cofold(char *string, char *constraints = NULL, float *OUTPUT);
%ignore cofold;
%include  "../H/cofold.h"

//%subsection "Partition function Folding"

%rename (pf_fold) my_pf_fold;
%{
  char *my_pf_fold(char *string, char *constraints, float *energy) {
    char *struc;
    struc = calloc(strlen(string)+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(string));
    *energy = pf_fold(string, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_pf_fold;
char *my_pf_fold(char *string, char *constraints = NULL, float *OUTPUT);

%ignore pf_fold;
%include  "../H/part_func.h"

%newobject pbacktrack;
extern char *pbacktrack(char *sequence);

//%subsection "Inverse Folding"

%rename (inverse_fold) my_inverse_fold;
%{
  char *my_inverse_fold(char *start, const char *target, float *cost) {
    char *seq;
    int n;
    n = strlen(target);
    seq = random_string(n, symbolset);
    if (start)
      strncpy(seq, start, strlen(start));
    *cost = inverse_fold(seq, target);
    if (start)
      /* for backward compatibility modify start */
      strncpy(start, seq, strlen(start));
    return(seq);
  }
%}

%newobject my_inverse_fold;
char * my_inverse_fold(char *start, const char *target, float *OUTPUT);

%rename (inverse_pf_fold) my_inverse_pf_fold;
%{
  char *my_inverse_pf_fold(char *start, const char *target, float *cost) {
    char *seq;
    int n;
    n = strlen(target);
    seq = random_string(n, symbolset);
    if (start) strncpy(seq, start, n);
    *cost = inverse_pf_fold(seq, target);
    if (start)
      /* for backward compatibility modify start */
      strncpy(start, seq, strlen(start));
    return(seq);
  }
%}

%newobject my_inverse_pf_fold;
char * my_inverse_pf_fold(char *start, const char *target, float *OUTPUT);

%ignore inverse_fold;
%ignore inverse_pf_fold;
%include  "../H/inverse.h"

//%subsection "Global Variables to Modify Folding"
//extern double *pr;  /*  base pairing prob. matrix */
%include  "../H/fold_vars.h"
%extend bondT {
	bondT *get(int i) {
	   return self+i;
	}
}

%ignore alifold;
%include "../H/alifold.h"
%rename (alifold) my_alifold;

%{
  char *my_alifold(char **strings, char *constraints, float *energy) {
    char *struc;
    struc = calloc(strlen(strings[0])+1,sizeof(char));
    if (constraints && fold_constrained)
      strncpy(struc, constraints, strlen(strings[0]));
    *energy = alifold(strings, struc);
    if (constraints)
      strncpy(constraints, struc, strlen(constraints));
    return(struc);
  }
%}

%newobject my_alifold;
char *my_alifold(char **strings, char *constraints = NULL, float *OUTPUT);

%newobject consensus;
%newobject consensus_mis;
char *consensus(const char **AS);
char *consens_mis(const char **AS);


/*#######################################*/
/* Get coordinates for xy plot           */
/*#######################################*/

%{
  typedef struct {
    float X; /* X coords */
    float Y; /* Y coords */
  } COORDINATE;

  COORDINATE *get_xy_coordinates(const char *structure){
    int i;
    short *table = make_pair_table(structure);
    short length = (short) strlen(structure);

    COORDINATE *coords = (COORDINATE *) space((length+1)*sizeof(COORDINATE));
    float *X = (float *) space((length+1)*sizeof(float));
    float *Y = (float *) space((length+1)*sizeof(float));

    if (rna_plot_type == 0)
      simple_xy_coordinates(table, X, Y);
    else
      naview_xy_coordinates(table, X, Y);

    for(i=0;i<=length;i++){
      coords[i].X = X[i];
      coords[i].Y = Y[i]; 
    }
    free(table);
    free(X);
    free(Y);
    return(coords);
  }
%}

typedef struct {
  float X; /* X coords */
  float Y; /* Y coords */
} COORDINATE;

COORDINATE *get_xy_coordinates(const char *structure);

%extend COORDINATE {
  COORDINATE *get(int i) {
    return self+i;
  }

}

/*#################################*/
/* END get coordinates for xy plot */
/*#################################*/


//%include  "../H/subopt.h"
// from subopt.h

typedef struct {
  float energy;                            /* energy of structure */
  char *structure;
} SOLUTION;

%newobject subopt;
extern  SOLUTION *subopt (char *seq, char *constraint, int delta, FILE *fp=NULL);

extern  int subopt_sorted;                       /* sort output by energy */
%extend SOLUTION {
	SOLUTION *get(int i) {
//	   static int size=-1;
//	   if (size<0) {
//	     SOLUTION *s;
//	     for (s=self; s->structure; s++);
//	     size= (int) (s-self);
//	   }
//	   if (i>=size) {
//	     warn("value out of range");
//	     return NULL;
//	   }
	   return self+i;
	}

	int size() {
	   SOLUTION *s;
	   for (s=self; s->structure; s++);
	   return (int)(s-self);
	}

	~SOLUTION() {
	   SOLUTION *s;
	   for (s=self; s->structure; s++) free(s->structure);
	   free(self);
	}
}
%{
double get_pr(int i, int j) {
  int ii;
  if (i>j) {ii=i; i=j; j=ii;}
  return pr[iindx[i]-j];
}
%}
double get_pr(int i, int j);
/* Get probability of pair i.j from the pr array */

//%section "Parsing and Comparing Structures"
// from RNAstruct.h

%newobject b2HIT;
%newobject b2C;
%newobject b2Shapiro;
%newobject add_root;
%newobject expand_Shapiro;
%newobject expand_Full;
%newobject unexpand_Full;
%newobject unweight;
char *b2HIT(char *structure);     // Full   -> HIT    [incl. root]
char *b2C(char *structure);       // Full   -> Coarse [incl. root]
char *b2Shapiro(char *structure); // Full -> weighted Shapiro [i.r]
char *add_root(char *);           // {Tree} -> ({Tree}R)
char *expand_Shapiro(char *coarse); // add S for stacks to coarse struct
char *expand_Full(char *structure); // Full   -> FFull
char *unexpand_Full(char *ffull);   // FFull  -> Full
char *unweight(char *wcoarse);   // remove weights from coarse struct
void   unexpand_aligned_F(char *align[2]);
void   parse_structure(char *structure); // make structure statistics
int    loop_size[1000];       // loop sizes of a structure
int    helix_size[1000];      // helix sizes of a structure
int    loop_degree[1000];     // loop degrees of a structure
int    loops;                 // n of loops and stacks
int    unpaired, pairs;       // n of unpaired digits and pairs

%include  "../H/treedist.h"
%include  "../H/stringdist.h"
%newobject Make_bp_profile;
%include  "../H/profiledist.h"
// from dist_vars.h
int   edit_backtrack;  /* set to 1 if you want backtracking */
char *aligned_line[2]; /* containes alignment after backtracking */
int  cost_matrix;      /* 0 usual costs (default), 1 Shapiro's costs */

//%section "Utilities"
%newobject space;
%newobject time_stamp;
%newobject random_string;
%newobject get_line;
%newobject pack_structure;
%newobject unpack_structure;
%newobject make_pair_table;

%include "../H/utils.h"

// from read_epars.c
extern void  read_parameter_file(char *fname);
/* read energy parameters from file */
extern void  write_parameter_file(char *fname);
/* write energy parameters to file */

// this doesn't work currently
%inline %{
void *deref_any(void **ptr, int index) {
   /* dereference arbitray pointer */
   return (void *) ptr[index];
}
%}

// from params.h

extern paramT *scale_parameters(void);
extern paramT *copy_parameters(void);
extern paramT *set_parameters(paramT *dest);
%{
char *get_aligned_line(int i) {
  i = i % 2;
  return aligned_line[i];
}
%}

char *get_aligned_line(int);



//%include ptr2array.i



%inline %{
  short *make_loop_index(const char *structure) {
  /* number each position by which loop it belongs to (positions start at 0) */
    int i,hx,l,nl;
    int length;
    short *stack;
    short *loop;
    length = strlen(structure);
    stack = (short *) space(sizeof(short)*(length+1));
    loop = (short *) space(sizeof(short)*(length+2));
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
	  nrerror("unbalanced brackets in make_pair_table");
	}
      }
    }
    free(stack);
    return loop;
  }
%}

%inline %{
float energy_of_move(const char *string, char *structure, int mi, int mj) {
  extern int energy_of_struct_pt(const char *string, short * ptable,
				 short *s, short *s1);
#define ILLEGAL 999.;
  int i,j,hx,l,nl;
  int length;
  short *stack, *table, *loop;
  short *S;
  int energy;

  if (mj<0) {
    if ((structure[-mi]!='(') || (structure[-mj]!=')'))
      return 1001;  /* illegal delete pair */
  } else
    if ((structure[mi]!='.') || (structure[mj]!='.'))
      return 1002;  /* illegal add pair */

  /* make the pair table and loop index l*/
  length = strlen(structure);
  stack = (short *) space(sizeof(short)*(length+1));
  loop  = (short *) space(sizeof(short)*(length+2));
  table = (short *) space(sizeof(short)*(length+2));
  table[0] = length;
  hx=l=nl=0;
  for (i=1; i<=length; i++) {
    if (structure[i-1] == '(') {
      nl++; l=nl;
      stack[hx++]=i;
    }
    loop[i]=l;
    if (structure[i-1] ==')') {
      j=stack[--hx];
      if (hx>0)
	l = loop[stack[hx-1]];  /* index of enclosing loop   */
      else l=0;                 /* external loop has index 0 */
      if (hx<0) {
	fprintf(stderr, "%s\n", structure);
	nrerror("unbalanced brackets in energy_of_move");
      }
      table[i]=j;
      table[j]=i;
    }
  }
  if (hx!=0) {
    fprintf(stderr, "%s\n", structure);
    nrerror("unbalanced brackets in energy_of_move");
  }

  if (loop[abs(mi)+1] != loop[abs(mj)+1]) { /* not in same loop => illegal */
    free(stack);
    free(loop);
    free(table);
    return 1003.;
  }

  /* if we get here the move is legal */
  if (mj<0) { /* delete pair */
    structure[-mi] = '.';
    structure[-mj] = '.';
    table[-mi+1] = table[-mj+1] = 0;
  } else { /* insert pair */
    structure[mi] = '(';
    structure[mj] = ')';
    table[mi+1] = mj+1;
    table[mj+1] = mi+1;
  }

  S = (short *) space(sizeof(short)*(length+1));
  S[0] = length;
  for (i=1; i<=length; i++) {
    char *pos;
    pos = strchr(Law_and_Order, string[i-1]);
    if (pos==NULL) S[i]=0;
    else S[i] = pos-Law_and_Order;
  }

  energy =  energy_of_struct_pt(string, table, S, S);

  free(S);
  free(stack);
  free(loop);
  free(table);
  return (float) energy/100.;
}
%}

%init %{
/* work around segfault when script tries to free symbolset */

symbolset = (char *) space(21);
strcpy(symbolset, "AUGC");

%}

//%include "../H/duplex.h"
typedef struct {
  int i;
  int j;
  char *structure;
  float energy;
} duplexT;

extern duplexT duplexfold(const char *s1, const char *s2);
extern duplexT aliduplexfold(const char **s1, const char **s2);
// problem: wrapper will not free the the structure hidden in duplexT
// even more problem for duplex_subopt()

%{
short *encode_seq(char *sequence) {
  unsigned int i,l;
  short *S;
  l = strlen(sequence);
  S = (short *) space(sizeof(short)*(l+2));
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

%include  "../H/PS_dot.h"
