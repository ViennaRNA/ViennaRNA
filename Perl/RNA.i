%module RNA
%{
#include  "../H/utils.h"
#include  "../H/fold_vars.h"
#undef fold
#include  "../H/fold.h"
#include  "../H/part_func.h"
#include  "../H/PS_dot.h"
#include  "../H/inverse.h"
#include  "../H/RNAstruct.h"
#include  "../H/treedist.h"
#include  "../H/stringdist.h"
#include  "../H/profiledist.h"
#include  "../H/dist_vars.h"
#include  "../H/pair_mat.h"
#include  "../H/subopt.h"
%}
//

%constant double VERSION = 0.3;
%include typemaps.i
%typemap(perl5,in) FILE * {
$target = IoIFP(sv_2io($source));
}

%title "Interface to the Vienna RNA library"
%section "Folding Routines"
%subsection "Minimum free Energy Folding"
%apply char *BOTH {char *structure};
%include  "../H/fold.h"
%subsection "Partition function Folding"
%include  "../H/part_func.h"
%clear char *structure;
%subsection "Inverse Folding"
%include  "../H/inverse.h"
%subsection "Global Variables to Modify Folding"
//extern double *pr;  /*  base pairing prob. matrix */
%include  "../H/fold_vars.h"
%include  "../H/subopt.h"

%addmethods SOLUTION {
	SOLUTION *get(int i) {
	   return self+i;
	}
	int size() {
	   SOLUTION *s;
	   for (s=self; s->structure; s++);
	   return (int)(s-self);
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

%section "Parsing and Comparing Structures"
// from RNAstruct.h

%new char *b2HIT(char *structure);   /* Full   -> HIT    [incl. root]      */
%new char *b2C(char *structure);     /* Full   -> Coarse [incl. root]      */
%new char *b2Shapiro(char *structure); /* Full -> weighted Shapiro [i.r.] */
%new char *add_root(char *);         /* {Tree} -> ({Tree}R)                */
%new char  *expand_Shapiro(char *coarse); // add S for stacks to coarse struct 
%new char  *expand_Full(char *structure); /* Full   -> FFull              */
%new char  *unexpand_Full(char *ffull);   /* FFull  -> Full               */
%new char  *unweight(char *wcoarse);   /* remove weights from coarse struct */
void   unexpand_aligned_F(char *align[2]);
void   parse_structure(char *structure); /* make structure statistics */
int    loop_size[1000];       /* loop sizes of a structure */
int    helix_size[1000];      /* helix sizes of a structure */
int    loop_degree[1000];     /* loop degrees of a structure */
int    loops;                  /* n of loops and stacks */
int    unpaired, pairs;        /* n of unpaired digits and pairs */

%include  "../H/treedist.h"
%include  "../H/stringdist.h"
%include  "../H/profiledist.h"
// from dist_vars.h
int   edit_backtrack;  /* set to 1 if you want backtracking */ 
char *aligned_line[2]; /* containes alignment after backtracking */
int  cost_matrix;      /* 0 usual costs (default), 1 Shapiro's costs */

%section "Utilities"
// from utils.h
%new void  *space(unsigned size);           /* allocate space safely */
void   nrerror(const char message[]);  /* die with error message */
void   init_rand(void);                /* make random number seeds */
short xsubi[3];               	       /* current 48bit random number */
double urn(void);                      /* random number from [0..1] */
int    int_urn(int from, int to);      /* random integer */
void   filecopy(FILE *from, FILE *to); /* inefficient `cp' */
%new char  *time_stamp(void);               /* current date in a string */
%new char  *random_string(int l, const char symbols[]);
/* random string of length l using characters from symbols[] */
int    hamming(const char *s1, const char *s2);
/* calculate hamming distance */
%new char  *get_line(FILE *fp); /* read one (arbitrary length) line from fp */
%new char *pack_structure(const char *struc);
/* pack secondary secondary structure, 5:1 compression using base 3 encoding */
%new char *unpack_structure(const char *packed);
/* unpack sec structure packed with pack_structure() */
%new short *make_pair_table(const char *structure);
/* returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
   0 if i is unpaired, table[0] contains the length of the structure. */
int bp_distance(const char *str1, const char *str2);
/* dist = {number of base pairs in one structure but not in the other} 
   same as edit distance with open-pair close-pair as move-set */

// from read_epars.c
extern void  read_parameter_file(char *fname);
/* read energy parameters from file */
extern void  write_parameter_file(char *fname);
/* write energy parameters to file */
//%include array.i
%include pointer.i
%inline %{
void *deref_any(void **ptr, int index) {
   /* dereference arbitray pointer */
   return (void *) ptr[index];
}
%}

%include ptr2array.i

%inline %{
  short *make_loop_index(const char *structure) {
  /* number each position by which loop it belongs to (positions start at 0) */
    int i,j,hx,l,nl;
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
  short *S, *S1;
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


// Convert Perl array reference int a char ** 
%typemap(perl5,in) char ** {
  AV *tempav;
  I32 len;
  int i;
  SV **tv;
  if (!SvROK($source)) croak("$source is not a reference.");
  if (SvTYPE(SvRV($source)) != SVt_PVAV) croak("$source is not an array.");
  tempav = (AV*)SvRV($source);
  len = av_len(tempav);
  $target = (char **) malloc((len+2)*sizeof(char *));
  for (i = 0; i <= len; i++) {
    tv = av_fetch(tempav, i, 0);
    $target[i] = (char *) SvPV(*tv,PL_na);
  }
  $target[i] = 0;
}
// This cleans up our char ** array after the function call
%typemap(perl5,freearg) char ** {
  free($source);
}
%include  "../H/PS_dot.h"

